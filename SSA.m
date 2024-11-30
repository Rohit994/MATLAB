clc; clear;

% Define constants
K2 = 1;  % Fixed value for k2
k5 = 20.16e6; % Initial fixed value for k5
O2 = 31.2;    % Convert µM to M
H_plus = 0.302; % Convert µM to M
h = 0.31;
sites = 6.7e-5;

% Parameter ranges
k1_values = linspace(50, 100, 41);  % k1 from 50 to 100
k2_values = linspace(100, 300, 41); % k2 from 100 to 300
k5_values = linspace(10e6, 30e6, 41); % k5 from 15e6 to 25e6

% Initialize results
X_k1 = zeros(length(k1_values), 1);
X_k2 = zeros(length(k2_values), 1);
X_k5 = zeros(length(k5_values), 1);

% Time span for the ODE solution
tspan = [0, 2];  % 2 hours
X0 = 388;        % Initial condition for X

% Vary k1 while keeping k2 and k5 fixed
for i = 1:length(k1_values)
    k1 = k1_values(i);
    
    % Solve the ODE
    [t, X] = ode45(@(t, X) ode_function(t, X, K2, k1, K2, k5, O2, H_plus, h, sites), tspan, X0);
    
    % Store the final concentration
    X_k1(i) = X(end);
end

% Vary k2 while keeping k1 and k5 fixed
for j = 1:length(k2_values)
    k2 = k2_values(j);
    
    % Solve the ODE
    [t, X] = ode45(@(t, X) ode_function(t, X, K2, k1_values(1), k2, k5, O2, H_plus, h, sites), tspan, X0);
    
    % Store the final concentration
    X_k2(j) = X(end);
end

% Vary k5 while keeping k1 and k2 fixed
for k = 1:length(k5_values)
    k5 = k5_values(k);
    
    % Solve the ODE
    [t, X] = ode45(@(t, X) ode_function(t, X, K2, k1_values(1), K2, k5, O2, H_plus, h, sites), tspan, X0);
    
    % Store the final concentration
    X_k5(k) = X(end);
end

% Plotting
figure;

% Plot for k1 effect
subplot(3, 1, 1);
plot(k1_values, X_k1, '-o', 'LineWidth', 1.5, 'Color', 'b');
xlabel('k1 (M^{-1}s^{-1})');
ylabel('X(40)');
title('Effect of k1 on Concentration at the End of 2 Hours');
grid on;

% Plot for k2 effect
subplot(3, 1, 2);
plot(k2_values, X_k2, '-o', 'LineWidth', 1.5, 'Color', 'g');
xlabel('k2 (M^{-1}s^{-1})');
ylabel('X(40)');
title('Effect of k2 on Concentration at the End of 2 Hours');
grid on;

% Plot for k5 effect
subplot(3, 1, 3);
plot(k5_values, X_k5, '-o', 'LineWidth', 1.5, 'Color', 'r');
xlabel('k5 (M^{-1}s^{-1})');
ylabel('X(40)');
title('Effect of k5 on Concentration at the End of 2 Hours');
grid on;

% Display y-axis in scientific notation for clarity
set(gca, 'YScale', 'linear');  % Set y-axis scale to linear
yticks = get(gca, 'YTick');     % Get current y-ticks
set(gca, 'YTickLabel', arrayfun(@(x) sprintf('%.2e', x), yticks, 'UniformOutput', false));  % Format y-ticks

% Change x-axis ticks for the third plot
subplot(3, 1, 3); % Go back to the third subplot
xticks(15e6:1e6:25e6); % Set ticks from 15e6 to 25e6

% Define the ODE function (constants are passed as parameters)
function dXdt = ode_function(t, ~, K2, k1, k2, k5, O2, H_plus, h, sites)
    % Term calculations for readability
    term1 = -(k2 * (-0.532 + (2.469 / (1 + exp((t - 1.891) / 0.299))))) * h;
    term2 = k1 * K2 * sites * O2 * (-0.532 + (2.469 / (1 + exp((t - 1.891) / 0.299)))) * H_plus;
    term3 = k5 * (0.0532 - 0.0553 * t + 3.088 * t^2 - 1.491 * t^3) * (-4.695 + 8.228e-4 * exp((t - 0.123) / 0.220));
    
    % Differential equation expression
    dXdt = (term1 + term2 + term3) * 1e-6;
end
