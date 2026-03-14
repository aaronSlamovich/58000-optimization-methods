
% plotting griewank function
% x1 = linspace(-200,200,5000);
% x2 = linspace(-200,200,5000);
% f_x = 1/4000 .* (x1.^2 + x2.^2') - cos(x1/sqrt(1))' * cos(x2/sqrt(2)) + 1;
%mesh(f_x)

% returns value of griewank function at (x1,x2)
function f_x = griewank(cp)
    f_x = 1/4000 .* (cp(1).^2 + cp(2).^2') - cos(cp(1)/sqrt(1))' * cos(cp(2)/sqrt(2)) + 1;
end

% returns value of derivative of griewank function at current point 'cp'
function [d_f_x] = d_griewank(cp)
    d_f_x(1) = 2/4000 .* cp(1) + 1/sqrt(1) * sin(cp(1)/sqrt(1))' * cos(cp(2)/sqrt(2)); % yes, I know 1/sqrt(1) = 1, I am just writing my code explicitly
    d_f_x(2) = 2/4000 .* cp(2) + 1/sqrt(2) * cos(cp(1)/sqrt(1))' * sin(cp(2)/sqrt(2));
end

function [a0,b0] = bracket(cp, d_cp)
    % find bracket for golden section
    %phi = @(alpha) griewank(cp - alpha*d_cp);
    a0 = 0;

    space = 0.075;
    
    x0 = cp;
    
    dx = d_cp;
    
    x1 = x0 - space * dx;
    space = 2 * space;
    x2 = x1 - space * dx;
    space = 2 * space;
    
    while griewank(x1) >= griewank(x0)
        x1 = x1 - space * dx;
        space = space / 2;
        x2 = x1;
    end
    
    while griewank(x2) <= griewank(x1)
        x2 = x2 - space * dx;
        space = 2 * space;
    end
    b0 = norm((x2-x0) ./ d_cp);
end

function alpha = golden_search(a0, b0, cp, d_cp)
    % a0,b0 are scalars returned by the bracketing
    phi = @(alpha) griewank(cp - alpha*d_cp);
    r = (3 - sqrt(5)) / 2; % r for ratio
    
    % remember log property log_b(x) = log_a(x) / log_a(b)
    n = ceil(log(0.23/(b0-a0)) / log(1-r)); % finds iterations n S.T. uncertainty < 0.23
    
    for i=1:n
        
        a1 = a0 + r * (b0 - a0);
        b1 = b0 - r * (b0 - a0);
    
        % choose the greater between f(a1) and f(b1)
        if (phi(a1)) > (phi(b1))
            a0 = a1;
        else
            b0 = b1;
        end
        
    end
    if phi(a1) < phi(b1)
        alpha = a1;
    else 
        alpha = b1;
    end
end

% ----------------steepest descent algorithm---------------
cp = [100, 150]; % cp stands for 'current point' point of SD alg

for i = 1:100
    % calculate derivative at point
    cp_val = griewank(cp); % evaluate griewank at cp
    d_cp = d_griewank(cp); % evaluate derivative of griewank at cp
    fprintf('griewank(%d,%d) = %d\n', cp(1), cp(2), cp_val)
    fprintf('derivative of griewank(%d,%d) = (%d,%d)\n', cp(1), cp(2), d_cp(1), d_cp(2))
    
    % find bracket for golden section
    [x1,x2] = bracket(cp, d_cp);
    
    % conduct golden section search in direction of negative gradient
    
    alpha = golden_search(x1, x2, cp, d_cp);
    cp = cp - alpha * d_cp;
end