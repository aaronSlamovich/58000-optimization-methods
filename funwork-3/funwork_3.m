clc;
clear;
close all;

start1 = [20, 30];
start2 = [-100, -150];

maxit = 10;

function f_x = griewank(cp)
    f_x = 1/4000 * (cp(1)^2 + cp(2)^2) - cos(cp(1)) * cos(cp(2)/sqrt(2)) + 1;
end

function d_f_x = d_griewank(cp)
    d_f_x(1) = cp(1)/2000 + sin(cp(1)) * cos(cp(2)/sqrt(2));
    d_f_x(2) = cp(2)/2000 + (1/sqrt(2)) * cos(cp(1)) * sin(cp(2)/sqrt(2));
end

function [a0,b0] = bracket(cp, d_cp)
    phi = @(alpha) griewank(cp - alpha*d_cp);

    a0 = 0;
    step = 10;
    a1 = step;

    while phi(a1) >= phi(a0) && step > 1e-12 
        step = step / 2;
        a1 = step;
    end

    if step <= 1e-12
        b0 = step;
        return;
    end

    a2 = 2*a1;
    while phi(a2) <= phi(a1)
        a1 = a2;
        a2 = 2*a2;
    end

    a0 = 0;
    b0 = a2;
    disp(a0)
    disp(b0)
end

function alpha = golden_search(a0, b0, cp, d_cp)
    phi = @(alpha) griewank(cp - alpha*d_cp);
    r = (3 - sqrt(5)) / 2;

    tol = 1e-4;
    n = ceil(log(tol/(b0-a0)) / log(1-r));
    n = max(n,1);

    for i = 1:n
        a1 = a0 + r * (b0 - a0);
        b1 = b0 - r * (b0 - a0);

        if phi(a1) > phi(b1)
            a0 = a1;
        else
            b0 = b1;
        end
    end

    alpha = (a0 + b0)/2;
    disp('alpha:')
    %disp(alpha)
end


function plot_level_sets(path1, path2, plot_name)
    x1 = linspace(-120,120,300);
    x2 = linspace(-120,120,300);
    [X1,X2] = meshgrid(x1,x2);
    Z = 1/4000 .* (X1.^2 + X2.^2) - cos(X1) .* cos(X2./sqrt(2)) + 1;

    figure;
    tiledlayout(1,2);

    % ---------- first starting point ----------
    nexttile;
    contour(X1, X2, Z, 35, 'LineWidth', 1);
    hold on;
    axis equal;
    grid on;
    title([plot_name '  (start 1)']);
    xlabel('x_1');
    ylabel('x_2');

    plot(path1(:,1), path1(:,2), 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    for k = 1:size(path1,1)-1
        quiver(path1(k,1), path1(k,2), ...
               path1(k+1,1)-path1(k,1), path1(k+1,2)-path1(k,2), ...
               0, 'r', 'LineWidth', 1.2, 'MaxHeadSize', 0.5);
    end
    scatter(path1(1,1), path1(1,2), 70, 'g', 'filled');
    scatter(path1(end,1), path1(end,2), 70, 'm', 'filled');
    hold off;

    % ---------- second starting point ----------
    nexttile;
    contour(X1, X2, Z, 35, 'LineWidth', 1);
    hold on;
    axis equal;
    grid on;
    title([plot_name '  (start 2)']);
    xlabel('x_1');
    ylabel('x_2');

    plot(path2(:,1), path2(:,2), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    for k = 1:size(path2,1)-1
        quiver(path2(k,1), path2(k,2), ...
               path2(k+1,1)-path2(k,1), path2(k+1,2)-path2(k,2), ...
               0, 'b', 'LineWidth', 1.2, 'MaxHeadSize', 0.5);
    end
    scatter(path2(1,1), path2(1,2), 70, 'g', 'filled');
    scatter(path2(end,1), path2(end,2), 70, 'm', 'filled');
    hold off;
end

% -----------------Steepest Descent----------------------

% ----- starting point 1 -----
cp = start1;
path_sd_1 = zeros(maxit+1,2);
path_sd_1(1,:) = cp;

for i = 1:maxit
    disp(griewank(cp))
    cp_val = griewank(cp);
    d_cp = d_griewank(cp);
    dir_cp = d_cp / norm(d_cp);

    fprintf('griewank(%.6f,%.6f) = %.6f\n', cp(1), cp(2), cp_val)
    fprintf('direction = (%.6f,%.6f)\n', d_cp(1), d_cp(2))

    [x1,x2] = bracket(cp, dir_cp);
    alpha = golden_search(x1, x2, cp, dir_cp);
    cp = cp - alpha * dir_cp;

    path_sd_1(i+1,:) = cp;
end

% ----- starting point 2 -----
cp = start2;
path_sd_2 = zeros(maxit+1,2);
path_sd_2(1,:) = cp;

for i = 1:maxit
    disp(griewank(cp))
    cp_val = griewank(cp);
    d_cp = d_griewank(cp);
    dir_cp = d_cp / norm(d_cp);

    fprintf('griewank(%.6f,%.6f) = %.6f\n', cp(1), cp(2), cp_val)
    fprintf('direction = (%.6f,%.6f)\n', d_cp(1), d_cp(2))

    [x1,x2] = bracket(cp, dir_cp);
    alpha = golden_search(x1, x2, cp, dir_cp);
    cp = cp - alpha * dir_cp;

    path_sd_2(i+1,:) = cp;
end

plot_level_sets(path_sd_1, path_sd_2, 'Steepest Descent');

% ----------------Conjugate Gradient-------------------

% ----- starting point 1 -----
cp = start1;
path_cg_1 = zeros(maxit+1,2);
path_cg_1(1,:) = cp;
g_prev = [];
p_prev = [];

for i = 1:maxit
    disp(griewank(cp))
    cp_val = griewank(cp);
    g = d_griewank(cp);

    if i == 1
        d_cp = g;
    else
        beta = (g*g') / (g_prev*g_prev');
        if mod(i-1, 2) == 0
            beta = 0;
        end
        d_cp = g + beta * p_prev;
    end

    fprintf('griewank(%.6f,%.6f) = %.6f\n', cp(1), cp(2), cp_val)
    fprintf('direction = (%.6f,%.6f)\n', d_cp(1), d_cp(2))

    [x1,x2] = bracket(cp, d_cp);
    alpha = golden_search(x1, x2, cp, d_cp);
    cp = cp - alpha * d_cp;

    g_prev = g;
    p_prev = d_cp;
    path_cg_1(i+1,:) = cp;
end

% ----- starting point 2 -----
cp = start2;
path_cg_2 = zeros(maxit+1,2);
path_cg_2(1,:) = cp;
g_prev = [];
p_prev = [];

for i = 1:maxit
    disp(griewank(cp))
    cp_val = griewank(cp);
    g = d_griewank(cp);

    if i == 1
        d_cp = g;
    else
        beta = (g*g') / (g_prev*g_prev');
        if mod(i-1, 2) == 0
            beta = 0;
        end
        d_cp = g + beta * p_prev;
    end

    fprintf('griewank(%.6f,%.6f) = %.6f\n', cp(1), cp(2), cp_val)
    fprintf('direction = (%.6f,%.6f)\n', d_cp(1), d_cp(2))

    [x1,x2] = bracket(cp, d_cp);
    alpha = golden_search(x1, x2, cp, d_cp);
    cp = cp - alpha * d_cp;

    g_prev = g;
    p_prev = d_cp;
    path_cg_2(i+1,:) = cp;
end

plot_level_sets(path_cg_1, path_cg_2, 'Powell Conjugate Gradient');

% -----------------------rank one----------------------

% ----- starting point 1 -----
cp = start1;
path_rank1_1 = zeros(maxit+1,2);
path_rank1_1(1,:) = cp;
H = eye(2);

for i = 1:maxit
    disp(griewank(cp))
    cp_val = griewank(cp);
    g = d_griewank(cp);
    d_cp = H * g';
    d_cp = d_cp';

    fprintf('griewank(%.6f,%.6f) = %.6f\n', cp(1), cp(2), cp_val)
    fprintf('direction = (%.6f,%.6f)\n', d_cp(1), d_cp(2))

    [x1,x2] = bracket(cp, d_cp);
    alpha = golden_search(x1, x2, cp, d_cp);

    cp_old = cp;
    cp_new = cp_old - alpha * d_cp;
    g_new = d_griewank(cp_new);

    s = (cp_new - cp_old)';
    y = (g_new - g)';
    v = s - H*y;

    if abs(v' * y) > 1e-12
        H = H + (v * v') / (v' * y);
    end

    cp = cp_new;
    path_rank1_1(i+1,:) = cp;
end

% ----- starting point 2 -----
cp = start2;
path_rank1_2 = zeros(maxit+1,2);
path_rank1_2(1,:) = cp;
H = eye(2);

for i = 1:maxit
    disp(griewank(cp))
    cp_val = griewank(cp);
    g = d_griewank(cp);
    d_cp = H * g';
    d_cp = d_cp';

    fprintf('griewank(%.6f,%.6f) = %.6f\n', cp(1), cp(2), cp_val)
    fprintf('direction = (%.6f,%.6f)\n', d_cp(1), d_cp(2))

    [x1,x2] = bracket(cp, d_cp);
    alpha = golden_search(x1, x2, cp, d_cp);

    cp_old = cp;
    cp_new = cp_old - alpha * d_cp;
    g_new = d_griewank(cp_new);

    s = (cp_new - cp_old)';
    y = (g_new - g)';
    v = s - H*y;

    if abs(v' * y) > 1e-12
        H = H + (v * v') / (v' * y);
    end

    cp = cp_new;
    path_rank1_2(i+1,:) = cp;
end

plot_level_sets(path_rank1_1, path_rank1_2, 'Rank One Correction');

%---------------DFP-----------------------

% ----- starting point 1 -----
cp = start1;
path_dfp_1 = zeros(maxit+1,2);
path_dfp_1(1,:) = cp;
H = eye(2);

for i = 1:maxit
    disp(griewank(cp))
    cp_val = griewank(cp);
    g = d_griewank(cp);
    d_cp = H * g';
    d_cp = d_cp';

    fprintf('griewank(%.6f,%.6f) = %.6f\n', cp(1), cp(2), cp_val)
    fprintf('direction = (%.6f,%.6f)\n', d_cp(1), d_cp(2))

    [x1,x2] = bracket(cp, d_cp);
    alpha = golden_search(x1, x2, cp, d_cp);

    cp_old = cp;
    cp_new = cp_old - alpha * d_cp;
    g_new = d_griewank(cp_new);

    s = (cp_new - cp_old)';
    y = (g_new - g)';

    if abs(s' * y) > 1e-12 && abs(y' * H * y) > 1e-12
        H = H + (s*s')/(s'*y) - (H*y*y'*H)/(y'*H*y);
    end

    cp = cp_new;
    path_dfp_1(i+1,:) = cp;
end

% ----- starting point 2 -----
cp = start2;
path_dfp_2 = zeros(maxit+1,2);
path_dfp_2(1,:) = cp;
H = eye(2);

for i = 1:maxit
    disp(griewank(cp))
    cp_val = griewank(cp);
    g = d_griewank(cp);
    d_cp = H * g';
    d_cp = d_cp';

    fprintf('griewank(%.6f,%.6f) = %.6f\n', cp(1), cp(2), cp_val)
    fprintf('direction = (%.6f,%.6f)\n', d_cp(1), d_cp(2))

    [x1,x2] = bracket(cp, d_cp);
    alpha = golden_search(x1, x2, cp, d_cp);

    cp_old = cp;
    cp_new = cp_old - alpha * d_cp;
    g_new = d_griewank(cp_new);

    s = (cp_new - cp_old)';
    y = (g_new - g)';

    if abs(s' * y) > 1e-12 && abs(y' * H * y) > 1e-12
        H = H + (s*s')/(s'*y) - (H*y*y'*H)/(y'*H*y);
    end

    cp = cp_new;
    path_dfp_2(i+1,:) = cp;
end

plot_level_sets(path_dfp_1, path_dfp_2, 'DFP');

% ----------------BFGS--------------------

% ----- starting point 1 -----
cp = start1;
path_bfgs_1 = zeros(maxit+1,2);
path_bfgs_1(1,:) = cp;
H = eye(2);

for i = 1:maxit
    disp(griewank(cp))
    cp_val = griewank(cp);
    g = d_griewank(cp);
    d_cp = H * g';
    d_cp = d_cp';

    fprintf('griewank(%.6f,%.6f) = %.6f\n', cp(1), cp(2), cp_val)
    fprintf('direction = (%.6f,%.6f)\n', d_cp(1), d_cp(2))

    [x1,x2] = bracket(cp, d_cp);
    alpha = golden_search(x1, x2, cp, d_cp);

    cp_old = cp;
    cp_new = cp_old - alpha * d_cp;
    g_new = d_griewank(cp_new);

    s = (cp_new - cp_old)';
    y = (g_new - g)';
    rho = y' * s;

    if abs(rho) > 1e-12
        I = eye(2);
        H = (I - (s*y')/rho) * H * (I - (y*s')/rho) + (s*s')/rho;
    end

    cp = cp_new;
    path_bfgs_1(i+1,:) = cp;
end

% ----- starting point 2 -----
cp = start2;
path_bfgs_2 = zeros(maxit+1,2);
path_bfgs_2(1,:) = cp;
H = eye(2);

for i = 1:maxit
    disp(griewank(cp))
    cp_val = griewank(cp);
    g = d_griewank(cp);
    d_cp = H * g';
    d_cp = d_cp';

    fprintf('griewank(%.6f,%.6f) = %.6f\n', cp(1), cp(2), cp_val)
    fprintf('direction = (%.6f,%.6f)\n', d_cp(1), d_cp(2))

    [x1,x2] = bracket(cp, d_cp);
    alpha = golden_search(x1, x2, cp, d_cp);

    cp_old = cp;
    cp_new = cp_old - alpha * d_cp;
    g_new = d_griewank(cp_new);

    s = (cp_new - cp_old)';
    y = (g_new - g)';
    rho = y' * s;

    if abs(rho) > 1e-12
        I = eye(2);
        H = (I - (s*y')/rho) * H * (I - (y*s')/rho) + (s*s')/rho;
    end

    cp = cp_new;
    path_bfgs_2(i+1,:) = cp;
end

plot_level_sets(path_bfgs_1, path_bfgs_2, 'BFGS');