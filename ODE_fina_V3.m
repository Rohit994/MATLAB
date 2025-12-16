clc; clear;
% Define constants outside the function

k0 = 6;
TiO2 = 1e-5;
I0 = 8.28e-7;
k1 = 183.6;
k2 = 97.2;
k5 = 36e6;   %µM-1h-1
O2 = 31.2;    % µM
H_plus = 0.03; %µM
h= 0.3; %µM


% Define the ODE function (constants are passed as parameters)
function dXdt = ode_function(t, ~, k0, TiO2, I0, k1, k2, k5, O2, H_plus, h)
    % Term calculations for readability
    CR = -4.695e-4 + 0.432e-8 * exp((t-0.123)./0.220);
    e_cb = -0.003 + (0.012./(1 + exp((t-1.891)/0.299)));
    OH = 2.659e-4 - 2.766e-4 * t + 0.015 * t^2 - 0.007 * t^3;
    dOHdt = -2.766e-4 + 0.03*t - 0.021*t^2;
   
    term1 = -k0 * I0 * TiO2;
    term2 = k2 * (e_cb) * h;
    term3 = -k1 * O2 * e_cb * H_plus;
    term4 = k5 * (CR)  * (OH);


    % Differential equation expression
    dXdt = (term1 + term2 + term3 + term4 + dOHdt);
end

% Initial condition and time span
X0 = 400;            % Initial condition for X
tspan = [0, 2];  % Time span for the ODE solution

% Solve the ODE using ode45, passing constants to the function
[t, X] = ode45(@(t, X) ode_function(t, X, k0, TiO2, I0, k1, k2, k5, O2, H_plus, h), tspan, X0);

time = [2 1.5 1.0 0.5 0.0];
Conc_1 = [173.34, 301.98, 359.60, 381.91, 388.31];
Conc_2 = [185.03, 301.17, 366.26, 371.68, 369.88];
Conc_3 = [179.19, 303.58, 358.94, 370.80, 378.10];
% Plot the results
figure(1);
plot(t, X, 'LineWidth', 2);
xlabel('Time (t)');
ylabel('x(t)');
title('Solution of dx/dt using ode45');
grid on;
hold
scatter(time, Conc_1, 'filled', 'MarkerEdgeColor', 'b'); % Conc.1 in blue
hold on;
scatter(time, Conc_2, 'filled', 'MarkerEdgeColor', 'r'); % Conc.2 in red
hold on;
scatter(time, Conc_3, 'filled', 'MarkerEdgeColor', 'r'); % Conc.2 in red
hold on;
