clear; clc; close all;

% derived value
analytical_val = 23 * pi / 192;

% MC set-up 
N_vals = round(logspace(2, 7, 100));
mc_estimates = zeros(size(N_vals));

% MC sim
for i = 1 : length(N_vals)
    N = N_vals(i);

    % make N rand vals for each X/Y/Z
    X = rand(N, 1);
    Y = rand(N, 1);
    Z = rand(N, 1);

    % bool calculation
    mc_bool = (X.^2 + Y.^2 < Z) & (Z.^2 > X .* Y);
    % average results
    mc_estimates(i) = sum(mc_bool) / N;
end

% plotting
fig = figure();
fig.Theme = 'light';
hold on;
semilogx(N_vals, mc_estimates, 'LineWidth', 1.5)
yline(analytical_val, 'r--', 'LineWidth', 2, 'Label', 'Analytical Solution')

grid on;
set(gca, 'XScale', 'log')
xlabel('Sample Size')
ylabel('Estimated Probability')
title('HW1: Monte Carlo Convergence')
legend('Monte Carlo Estimate', 'Analytical Solution', 'Location', 'best');

print(fig, 'APPM4560HW1_SimFig', '-dpng', '-r300')