% ---------minimizing griewank function (Problem 2)----------------

% returns vector of solutions for a matrix of nx2 coordinates
function output = griewank(x)
    output = 1 + x(:,1).^2 ./ 4000 + x(:,2).^2 ./ 4000 - cos(x(:,1)).*cos(x(:,2)/sqrt(2));
end

% returns best overall corrdinate and vector of best coordinates for each
% particle
function [gbest, pbest] = evaluate(swarm, pbest)
    swarm_outputs = griewank(swarm);
    pbest_outputs = griewank(pbest);
    pbest(swarm_outputs < pbest_outputs, :) = swarm(swarm_outputs < pbest_outputs, :);
    pbest_outputs = griewank(pbest);
    [mini, min_idx] = min(pbest_outputs);
    gbest = pbest(min_idx,:);
end

num_particles = 5;

% make a 5x5x2 matrix first page = x1 second page = x2
swarm = rand([num_particles 2]); % generate 25 random particles
swarm = swarm - 0.5; % allow results to be negative
swarm = swarm * 14; % scale up to -7x7

[gbest, pbest] = evaluate(swarm, swarm);
w = 0.9;
c1 = 2;
c2 = 2;
v_old = [0,0];
best = [];
average = [];
worst = [];

for i = 1:100
    % generate random vectors
    r = rand([num_particles, 2]);
    s = rand([num_particles, 2]);

    % calculate new velocity
    v_new = w .* v_old + c1 .* r .* (pbest - swarm) + c2 .* s .* (gbest - swarm);

    %update swarm
    swarm = swarm + v_new;
    v_old = v_new;

    [gbest, pbest] = evaluate(swarm, pbest);

    % for plotting
    outputs = griewank(swarm);
    best = [best griewank(gbest)];
    average = [average mean(griewank(pbest))];
    worst = [worst max(griewank(pbest))];
end

% plotting
% figure;
% plot(best, 'LineWidth', 2); hold on;
% plot(average, 'LineWidth', 2);
% plot(worst, 'LineWidth', 2);
% grid on;
% xlabel('Iteration');
% ylabel('Objective value');
% title('PSO Performance on Griewank Function');
% legend('Global Best', 'Average', 'Worst', 'Location', 'best');
% hold off;


% ---------minimizing peaks function (Problem 3)---------------

function output = my_peaks(x)
    a = 3.*(1-x(:,1)).^2.*exp(-1*x(:,1).^2-(x(:,2)+1).^2);
    b = 10*(x(:,1)./5 - x(:,1).^3 - x(:,2).^5).*exp(-x(:,1).^2-x(:,2).^2);
    c = (exp(-1*(x(:,1)+1).^2 - x(:,2).^2)) ./3;
    output = a-b-c;
end

% returns best overall corrdinate and vector of best coordinates for each
% particle
function [gbest, pbest] = evaluate_peaks(swarm, pbest)
    swarm_outputs = my_peaks(swarm);
    pbest_outputs = my_peaks(pbest);
    pbest(swarm_outputs < pbest_outputs, :) = swarm(swarm_outputs < pbest_outputs, :);
    pbest_outputs = my_peaks(pbest);
    [mini, min_idx] = min(pbest_outputs);
    gbest = pbest(min_idx,:); % change to preserve current gbest if new generation doesn't improve
end

num_particles = 5;

% make a 5x5x2 matrix first page = x1 second page = x2
swarm = rand([num_particles 2]); % generate 25 random particles
swarm = swarm - 0.5; % allow results to be negative
swarm = swarm * 6; % scale up to -3x3

[gbest, pbest] = evaluate(swarm, swarm);
w = 0.9;
c1 = 2;
c2 = 2;
v_old = [0,0];
best = [];
average = [];
worst = [];

for i = 1:100
    % generate random vectors
    r = rand([num_particles, 2]);
    s = rand([num_particles, 2]);

    % calculate new velocity
    v_new = w .* v_old + c1 .* r .* (pbest - swarm) + c2 .* s .* (gbest - swarm);

    %update swarm
    swarm = swarm + v_new;
    v_old = v_new;

    [gbest, pbest] = evaluate(swarm, pbest);

    % for plotting
    outputs = my_peaks(swarm);
    best = [best my_peaks(gbest)];
    average = [average mean(my_peaks(pbest))];
    worst = [worst max(my_peaks(pbest))];
end

% plotting
% figure;
% plot(best, 'LineWidth', 2); hold on;
% plot(average, 'LineWidth', 2);
% plot(worst, 'LineWidth', 2);
% grid on;
% xlabel('Iteration');
% ylabel('Objective value');
% title('PSO Performance on Peaks Function');
% legend('Global Best', 'Average', 'Worst', 'Location', 'best');
% hold off;

%------traveling salesperson genetic algorithm (Problem 4)-----------

% define city locations
x = [0.4306, 3.7094, 6.9330, 9.3582, 4.7758, 1.2910, 4.8381, 9.4560, 3.6774, 3.2849];
y = [7.7288, 2.9727, 1.7785, 6.9080, 2.6394, 4.5774, 8.43692, 8.8150, 7.0002, 7.5569];

loc = [x;y];

% initialize population
pop_size = 20;
P = [];
for i = 1:pop_size
    P = [P;randperm(10,10)];
end

% evaluate quality

% find distance for each chromosome
d = sum(abs(P(:,1:end-1) - P (:,2:end)),2);

% map distance to quality
q = 1./d; % might want to try a different mapping this might not be too good
q = q./sum(q); % probability of selection

% select mating pool with replacement
