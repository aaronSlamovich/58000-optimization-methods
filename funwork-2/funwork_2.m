% problem 2
disp('---------------problem 2 start------------------')

x = linspace(1,2,100);
f_x = 8 * exp(1-x) + 7 * log(x);

plot(x, f_x);
xlabel('x');
ylabel('f(x)');
title('Plot of f(x) = 8 * exp{1-x} + 7 * log(x)');

disp('is the function unimodal? Yes - we can see that it is monotonically decreasing for x =< 1.6 and monotonically increasing for x >= 1.6 within the bounds of our window')

disp('------------using golden section method-------------')
a0 = 1;
b0 = 2;

r = (3 - sqrt(5)) / 2; % r for ratio

% log_b(x) = log_a(x) / log_a(b)
n = ceil(log(0.23) / log(1-r));

for i=0:n
    
    a1 = a0 + r * (b0 - a0);
    b1 = b0 - r * (b0 - a0);

    fprintf('-----iteration =  %d  -----\n', i)
    fprintf('a_%d = %d\n', i, a0)
    fprintf('b_%d = %d\n', i, b0)
    fprintf('f(a_%d) = %d\n', i, 8 * exp(1-a0) + 7 * log(a0))
    fprintf('f(b_%d) = %d\n', i, 8 * exp(1-b0) + 7 * log(b0))
    fprintf('uncertainty interval = %d\n', b0-a0)

    % choose the greater between f(a1) and f(b1)
    if (8 * exp(1-a1) + 7 * log(a1)) > (8 * exp(1-b1) + 7 * log(b1))
        a0 = a1;
    else
        b0 = b1;
    end
    
end
    
disp('-------------------using fibonacci method-----------------')

a0 = 1;
b0 = 2;

% (1 + 2*0.05) / Fn+1 <= 0.23, N=10 is plenty

% log_b(x) = log_a(x) / log_a(b)
N = 5;

for i=0:N

    r = 1 - fibonacci(N-i+2) / fibonacci(N-i+3);
    disp(r)
    
    a1 = a0 + r * (b0 - a0);
    b1 = b0 - r * (b0 - a0);

    fprintf('-----iteration =  %d  -----\n', i)
    fprintf('a_%d = %d\n', i, a0)
    fprintf('b_%d = %d\n', i, b0)
    fprintf('f(a_%d) = %d\n', i, 8 * exp(1-a0) + 7 * log(a0))
    fprintf('f(b_%d) = %d\n', i, 8 * exp(1-b0) + 7 * log(b0))
    fprintf('uncertainty interval = %d\n', (b0-a0))

    % choose the greater between f(a1) and f(b1)
    if (8 * exp(1-a1) + 7 * log(a1)) > (8 * exp(1-b1) + 7 * log(b1))
        a0 = a1;
    else
        b0 = b1;
    end
    
end
   

% problem 3
disp('-------------------problem 3 start------------------')

% x = linspace(-0.001,0.005,1000);
% g_x = (2*x-1).^2 + 4*(4-1024*x).^4;

%plot(x, g_x)

x0 = 0;
x1 = 1;

while abs(x1-x0) >= abs(x0)*10e-5
    x2 = (g(x1)*x0 - g(x0)*x1) / (g(x1) - g(x0));
    x0 = x1;
    x1 = x2;
    disp(x2)
end

function g_x = g(x)
    g_x = (2*x-1).^2 + 4*(4-1024*x).^4;
end

disp('---------------problem 4 start-------------')

function f_x = f(x)
    f_x = x(1)^2 + x(1)*x(2) + x(2)^2;
end

function [dx1,dx2] = df(x)
    dx1 = 2 * x(1) + x(2);
    dx2 = 2*x(2) + x(1);
end

% bracketing
disp('------bracketing------')
space = 0.075;

x0 = [0.8, -0.25]';

dx = df(x0);
dx = dx / norm(dx);

x1 = x0 - space * dx;
space = 2 * space;
x2 = x1 - space * dx;
space = 2 * space;

while f(x1) >= f(x0) || f(x2) <= f(x1)
    x2 = x2 - space * dx;
    space = 2 * space;
end

disp('window start:')
disp(x0)
disp('window end:')
disp(x2)

% golden section
disp('-----golden section------')
a0 = x0;
b0 = x2;

r = (3 - sqrt(5)) / 2; % r for ratio

% log_b(x) = log_a(x) / log_a(b)
n = ceil(log(0.1) / log(1-r)) * (1/(norm(b0-a0)));

for i=0:n
    
    a1 = a0 + r * (b0 - a0);
    b1 = b0 - r * (b0 - a0);

    fprintf('-----iteration =  %d  -----\n', i)
    fprintf('f(a_%d) = %d\n', i, f(a0))
    fprintf('f(b_%d) = %d\n', i, f(b0))
    fprintf('uncertainty interval = %d\n', abs(b0-a0))

    % choose the greater between f(a1) and f(b1)
    if (f(a1) > f(b1))
        a0 = a1;
    else
        b0 = b1;
    end
    
end

disp('--------fibonacci--------')

a0 = x0;
b0 = x2;

%r = (3 - sqrt(5)) / 2; % r for ratio

% log_b(x) = log_a(x) / log_a(b)
n = 5;

for i=0:n
    r = 1 - fibonacci(n-i+2) / fibonacci(n-i+3);
    disp(r)
    a1 = a0 + r * (b0 - a0);
    b1 = b0 - r * (b0 - a0);

    fprintf('-----iteration =  %d  -----\n', i)
    fprintf('f(a_%d) = %d\n', i, f(a0))
    fprintf('f(b_%d) = %d\n', i, f(b0))
    fprintf('uncertainty interval = %d\n', abs(b0-a0))

    % choose the greater between f(a1) and f(b1)
    if (f(a1) > f(b1))
        a0 = a1;
    else
        b0 = b1;
    end
    
end

disp('-----------problem 5----------')

x = linspace(-2,2,100);
y = linspace(-1,3,100);

[X,Y] = meshgrid(x,y);

F = 100*(Y-X.^2).^2 + (1-X).^2;
figure;
surf(X,Y,F)
figure;
contour(X,Y,F)