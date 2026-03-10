
% returns value of griewank function at (x1,x2)
function f_x = griewank(cp)
    f_x = 1/4000 .* (cp(1).^2 + cp(2).^2') - cos(cp(1)/sqrt(1))' * cos(cp(2)/sqrt(2)) + 1;
end

% returns value of derivative of griewank function at current point 'cp'
function [d_f_x] = d_griewank(cp)
    d_f_x(1) = 2/4000 .* cp(1) - 1/sqrt(1) * sin(cp(1)/sqrt(1))' * cos(cp(2)/sqrt(2)); % yes, I know 1/sqrt(1) = 1, I am just writing my code explicitly
    d_f_x(2) = 2/4000 .* cp(2) - 1/sqrt(2) * cos(cp(1)/sqrt(1))' * sin(cp(2)/sqrt(2));
end

% plotting griewank function
% x1 = linspace(-200,200,5000);
% x2 = linspace(-200,200,5000);
% f_x = 1/4000 .* (x1.^2 + x2.^2') - cos(x1/sqrt(1))' * cos(x2/sqrt(2)) + 1;
%mesh(f_x)

% ----------------steepest descent algorithm---------------
cp = [100, 150]; % cp stands for 'current point' point of SD alg

% calculate derivative at point
cp_val = griewank(cp); % evaluate griewank at cp
d_cp = d_griewank(cp); % evaluate derivative of griewank at cp
fprintf('griewank(%d,%d) = %d\n', cp(1), cp(2), cp_val)
fprintf('derivative of griewank(%d,%d) = (%d,%d)\n', cp(1), cp(2), d_cp(1), d_cp(2))

% conduct line search in direction of negative gradient

