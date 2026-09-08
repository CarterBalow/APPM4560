clear; clc; close all;
%% PART B
num_samples = 10^5;

u = rand(1, num_samples);
z = 1 - sqrt(1- u);
z_mean = mean(z);

%% PART C

N = 10^4;
R = 10^3;
C = ones(1, R);

for n = 3 : (N - 1)
    C = C + (rand(1, R) < C / n);
end

z_sim = C / N;
emp_mean = mean(z_sim);
emp_std = std(z_sim);
emp_ratio = emp_std / emp_mean;

min_core = min(C);
max_core = max(C);

fig = figure();
fig.Theme = 'light';
hold on;

histogram(z_sim, 40, 'Normalization', 'pdf')

z_ax = linspace(0, 1, 500);
h_ax = 2 * (1 - z_ax);
plot(z_ax, h_ax, 'r-', 'LineWidth', 2)

xlim([0, 1]);
xlabel('z');
ylabel('Probability Density');
title('Core Node Ratio z = C/N Distribution');
subtitle(sprintf('N = %s, R = %s', num2str(N), num2str(R)));
legend('Empirical Probability (Simulation)', 'Theoretical h(z) = 2(1 - z)', 'Location', 'best');
grid on;

print(fig, 'IsotropicRedirectionGrowth', '-dpng', '-r300')