clear; clc; close all;

a = 0.9;
lamS = 10;
lamF = 1000;
N = 10^5;
t50 = 50E-3;

bool_fast = rand(N, 1) < a;
U = rand(N, 1);

T = zeros(N, 1);
T(bool_fast) = -log(U(bool_fast)) / lamF;
T(~bool_fast) = -log(U(~bool_fast)) / lamS;

theoretical_mean = (a / lamF) + ((1 - a) / lamS);
theoretical_p50 = a * exp(-lamF * t50) + (1 - a) * exp(-lamS * t50);

empirical_mean = mean(T);
empirical_p50 = mean(T > t50);

fig = figure();
fig.Theme = 'light';

hold on;

edges = linspace(0, max(T), 150); 
[counts, bin_edges] = histcounts(T, edges, 'Normalization', 'pdf');
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
bar(bin_centers, counts, 1);

t_grid = linspace(0, max(T), 1000);
f_theory = a * lamF * exp(-lamF * t_grid) + (1 - a) * lamS * exp(-lamS * t_grid);
plot(t_grid, f_theory, 'r-', 'LineWidth', 2);

set(gca, 'YScale', 'log');
ylim([1e-3, 1.2 * lamF]);
xlim([0, 0.6]);

xlabel('Dwell Time T (seconds)');
ylabel('Probability Density (log scale)');
title('Dwell Time Distribution');
legend('Empirical Dwell Time Probability', 'Theoretical Dwell Time Probability', 'Location', 'best');
grid on;

print(fig, 'DwellTimePDF', '-dpng', '-r300')