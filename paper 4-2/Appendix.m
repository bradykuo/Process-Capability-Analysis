function lcb(n,usl,r)
    % Initialize variables
    t = 0; b = 0; z = 0; e = 0; f = 0; g = 0;
    rt = zeros(2,1); x = 0; x1 = 0; x2 = 0; y = 0;
    y1 = 0; delta = 0; delta1 = 0; x3 = 0;
    
    % Data directly from Table 5
    data = [428.63 408.09 417.62 317.20 519.40 438.75 380.42 457.60 474.87 395.22 ...
            344.88 497.17 373.81 337.62 430.78 363.31 475.43 270.65 357.53 454.24 ...
            402.19 462.12 403.40 463.31 513.74 433.62 361.49 360.67 400.46 540.87 ...
            361.72 303.60 377.56 376.22 379.27 332.91 258.16 357.04 352.24 432.00 ...
            432.66 395.20 371.64 315.74 330.47 405.25 324.87 337.89 435.80 399.28 ...
            433.73 358.42 422.63 419.25 387.12 277.72 464.39 388.89 384.68 455.70 ...
            387.57 337.50 444.43 494.75 472.60 369.92 393.09 492.78 429.29 403.41 ...
            274.69 324.70 528.84 443.18 331.91 428.18 423.13 397.66 440.76 332.47]';

    % Compute X̄, S, bn-1, and C̃PU
    mdata = mean(data);
    stddata = std(data);
    
    % Calculate bn-1 using logarithms
    log_bn1 = 0.5 * log(2/(n-1)) + gammaln((n-1)/2) - gammaln((n-2)/2);
    b = exp(log_bn1);
    
    ecpu = (usl - mdata)/(3 * b * stddata);
    fprintf('The sample Mean is %g.\n',mdata)
    fprintf('The sample standard Deviation is %g.\n',stddata)
    fprintf('The Correction Factor is %g.\n',b)
    fprintf('The UMVUE of Cpu is %g.\n',ecpu)%%

    % Compute a good initial value of δ
    t = 3 * (ecpu/b) * sqrt(n);
    
    % Calculate R using logarithms
    log_R = 0.5 * log(2/(n-1)) + gammaln(n/2) - gammaln((n-1)/2);
    R = exp(log_R);
    
    z = norminv(r,0,1);
    e = R^2 - z^2 * (1 - R^2);
    f = -2 * R * (t * R^2 - z^2 * (1 - R^2) * t);
    g = (t * R^2 - z^2 * (1 - R^2) * t)^2 - z^2 * R^2 + z^4 * (1 - R^2);
    pol = [e f g];
    rt = roots(pol);
    x = nctinv(r, n-1, rt(2));
    x1 = nctinv(r, n-1, rt(2) + 0.01);
    x2 = (x1 - x)/9.5;
    y = abs(x - sqrt(n) * 3 * ecpu/b);
    y1 = y/x2;
    delta1 = rt(2) + y1 * 0.001;

    % Evaluate the lower confidence bound CU through numerical iterations
    delta = delta1 + 0.001 : 0.001 : delta1 + 0.5;
    for i = 1:1:500
        x3 = nctinv(r, n-1, delta(i));
        if(abs(x3 - sqrt(n) * 3 * ecpu/b)) <= 0.01
            cu = delta(i)/(3 * sqrt(n));
            fprintf('The Lower Confidence Bound is %g.\n',cu)
            break
        end
    end

    % Output the conclusive message
    fprintf('The true value of the process capability Cpu\n')
    fprintf('is no less than %g', cu)
    fprintf(' with %g',r)
    fprintf(' level of confidence\n')
end

% Call the function with the specified parameters
lcb(80, 650, 0.95);