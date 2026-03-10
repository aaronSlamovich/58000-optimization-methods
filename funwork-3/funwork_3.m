% plotting griewank function
x1 = linspace(-200,200,5000);
x2 = linspace(-200,200,5000);
f_x = 1/4000 .* (x1.^2 + x2.^2') - cos(x1/sqrt(1))' * cos(x2/sqrt(2)) + 1;
mesh(f_x)