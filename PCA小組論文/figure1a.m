% Define parameters
xi_values = 0:0.05:3.0; % ξ values from 0 to 3.0
beta = 0.05;
alpha = 0.05;
CLTPD = 1.00;
CAQL = 1.33; 

% Create meshgrid for plot (swap xi and C0)
C0_mesh = 1.166:0.0001:1.170;
xi_mesh = 0:0.05:3.0;  % Finer resolution for smoother surface
[C0, XI] = meshgrid(C0_mesh, xi_mesh);
N = zeros(size(C0));

% Calculate n values for each point in the mesh
for i = 1:size(C0,1)
   for j = 1:size(C0,2)
       xi = XI(i,j);
       
       % Calculate b1 and b2
       b1 = 3 * CAQL * (1 + xi^2)^0.5;
       b2 = 3 * CLTPD * (1 + xi^2)^0.5;
       
       % Define S1 and S2
       S1 = @(n, C0) integral(@(t) chi2cdf((b1.^2.*n)./(9.*C0.^2) - t.^2, n-1) .* ...
           (normpdf(t + xi.*sqrt(n)) + normpdf(t - xi.*sqrt(n))), ...
           0, (b1.*sqrt(n))./(3.*C0)) - (1-alpha);
       
       S2 = @(n, C0) integral(@(t) chi2cdf((b2.^2.*n)./(9.*C0.^2) - t.^2, n-1) .* ...
           (normpdf(t + xi.*sqrt(n)) + normpdf(t - xi.*sqrt(n))), ...
           0, (b2.*sqrt(n))./(3.*C0)) - beta;
       
       % Try to find n for this combination of xi and C0
       try
           obj_fun = @(n) [S1(n, C0(i,j)); S2(n, C0(i,j))];
           n_sol = fsolve(@(n) obj_fun(n), 50, optimset('Display', 'off'));
           if n_sol > 0 && n_sol < 200
               N(i,j) = n_sol;
           else
               N(i,j) = NaN;
           end
       catch
           N(i,j) = NaN;
       end
   end
end

% Create the surface plot with swapped axes
figure;
surf(C0, XI, N, 'EdgeColor', 'none');
colormap('jet');
colorbar;
shading interp;  % Smoother color transitions

% Customize the plot
xlabel('C_0');
ylabel('\xi');
zlabel('n');
grid on;
view(45, 35);  % Adjusted view angle

% Set axis limits
xlim([1.166 1.170]);
ylim([0 3.0]);
zlim([0 70]);

% Set specific axis ticks
set(gca, 'XDir', 'reverse');  % Reverse x-axis direction
set(gca, 'YTick', 0:0.5:3.0);  % Set y-axis ticks at 0.5 intervals