clc; clear;

% Define constants outside the function
K2 = 1;
k1 = 97.2;
k2 = 183.6;
k5 = 20.16e6;
O2 = 31.2;    % Convert µM to M
H_plus = 0.302; % Convert µM to M
h = 0.31;
sites = 6.7e-5;

% Define the range for v
v_values = linspace(68.1, 272.4, 5); % 5 evenly spaced values between 1e-6 and 5e-6

% Define the ODE function (constants are passed as parameters)
function dXdt = ode_function(t, ~, K2, k1, k2, k5, O2, H_plus, h, sites, v)
    % Term calculations for readability
    term1 = -(k2 * (-0.532 + (2.469 / (1 + exp((t - 1.891) / 0.299))))) * h;
    term2 = k1 * K2 * sites * O2 * (-0.532 + (2.469 / (1 + exp((t - 1.891) / 0.299)))) * H_plus;
    term3 = k5 * (0.0532 - 0.0553 * t + 3.088 * t^2 - 1.491 * t^3) * (-4.695 + 8.228e-4 * exp((t - 0.123) / 0.220));
    
    % Differential equation expression
    dXdt = 50e-6.*(term1 + term2 + term3)./v./2./0.35;
end

% Initial condition and time span
X0 = 400;            % Initial condition for X
tspan = [0, 2];      % Time span for the ODE solution

% Initialize figure
figure;
hold on;
colors = lines(length(v_values)); % Get a color map for the lines

% Loop over v values
for i = 1:length(v_values)
    v = v_values(i);
    
    % Solve the ODE using ode45
    [t, X] = ode45(@(t, X) ode_function(t, X, K2, k1, k2, k5, O2, H_plus, h, sites, v), tspan, X0);
    X
    % Plot the results
    plot(t, X, 'DisplayName', sprintf('v = %.1f l/h', v), 'Color', colors(i, :));
end

% Customize plot
xlabel('Time (h)');
ylabel('Concentration of X (μM)');
title('Concentration of X over Time for Different v Values');
legend show;
grid on;
%axis([0 2 0 400]); % [xmin xmax ymin ymax]
hold off;
