
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
    d_f_x(1) = 2/4000 .* cp(1) - 1/sqrt(1) * sin(cp(1)/sqrt(1))' * cos(cp(2)/sqrt(2)); % yes, I know 1/sqrt(1) = 1, I am just writing my code explicitly
    d_f_x(2) = 2/4000 .* cp(2) - 1/sqrt(2) * cos(cp(1)/sqrt(1))' * sin(cp(2)/sqrt(2));
end

function [x0,x2] = bracket(cp)
    % find bracket for golden section
    space = 0.075;
    
    x0 = cp;
    
    dx = d_cp;
    dx = dx / norm(dx);
    
    x1 = x0 - space * dx;
    space = 2 * space;
    x2 = x1 - space * dx;
    space = 2 * space;
    
    while griewank(x1) >= griewank(x0)
        x1 = x1 - space * dx;
        x2 = x1;
        space = 2 * space;
    end
    
    while griewank(x2) <= griewank(x1)
        x2 = x2 - space * dx;
        space = 2 * space;
    end
end

function new_point = golden_search(a0, b0)
    r = (3 - sqrt(5)) / 2; % r for ratio
    
    % remember log property log_b(x) = log_a(x) / log_a(b)
    n = ceil(log(0.23) / log(1-r)); % finds iterations n S.T. uncertainty < 0.23
    
    for i=0:n
        
        a1 = a0 + r * (b0 - a0);
        b1 = b0 - r * (b0 - a0);
    
        fprintf('-----iteration =  %d  -----\n', i)
        fprintf('a_%d = %d\n', i, a0)
        fprintf('b_%d = %d\n', i, b0)
        fprintf('f(a_%d) = %d\n', i, 8 * exp(1-a0) + 7 * log(a0))
        fprintf('f(b_%d) = %d\n', i, 8 * exp(1-b0) + 7 * log(b0))
        fprintf('uncertainty interval = %d\n', b0-a0)
    
        % this code assume a function is convex, which might not be a safe
        % assumption
        % choose the greater between f(a1) and f(b1)
        if (griewank(a1)) > (griewank(b1))
            a0 = a1;
        else
            b0 = b1;
        end
        
    end
    new_point = min(griewank(a1), griewank(b1));
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
    [x1,x2] = bracket(cp);
    
    % conduct golden section search in direction of negative gradient
    
    cp = golden_search(x1, x2);
end