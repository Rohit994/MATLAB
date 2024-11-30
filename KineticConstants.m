% Parameters
h = 0.309;         % Given value of h
%k3_nominal = 0.305;       % Nominal value of k3
k4_nominal = -0.81;      % Nominal value of k4

% Vary k3 and k4 by +/- 100%
k3_values = linspace(0.1, 0.5, 41);
k4_values = k4_nominal; % 50%, 100%, and 200% of nominal k4

% Initial condition
x0 = 400; % Given initial condition

% Time span (in hours)
tspan = [0, 2];  % Time interval from 0 to 2 hours


% Solve and plot for each combination of k3 and k4
figure;
hold on;
for k3 = k3_values
    for k4 = k4_values
        % Define the ODE with the current k3 and k4
        odefun = @(t, x) -x .* (-k3 * h - k4 * (0.053 - 0.553 * t + 3.088 * t^2 - 1.491 * t^3));
        
        % Solve the ODE
        [t, x] = ode45(odefun, tspan, x0);
        xres = x(41,:)
        % Plot the results
        plot(t, x, 'LineWidth', 1.5, 'DisplayName', sprintf('k3=%.3f, k4=%.2e', k3, k4));
        
    end
end
hold off;

% Add plot details
xlabel('Time t (hours)');
ylabel('Solution x(t)');
title('Solution of the ODE with Variations in k3 and k4');
legend('Location', 'best');
grid on;
