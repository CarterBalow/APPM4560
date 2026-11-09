clear; clc; close all;

Q = [-2 2 0 0;
    1 -3 2 0;
    0 2 -4 2;
    0 0 3 -3];
pi_ = [3 6 6 4] / 19;
x0 = [1 0 0 0];

lams = sort(real(eig(Q)), 'descend');
lam2 = lams(2);

xfunc = @(t) x0 * expm(t * Q);
tvfunc = @(t) 0.5 * sum(abs(xfunc(t) - pi_));

t = linspace(0, 5, 501);
X = zeros(numel(t), 4);
TV = zeros(numel(t), 1);
for k = 1 : numel(t)
    X(k, :) = xfunc(t(k));
    TV(k) = tvfunc(t(k));
end

tstar = fzero(@(s) tvfunc(s) - 0.01, [0.01 5]);
tpred = log(100) / abs(lam2);

C = TV(end) / exp(lam2 * t(end));

fig = figure();
fig.Theme = 'light';

subplot(1, 2, 1);
hold on;
cols = get(gca, 'ColorOrder');
for i = 1 : 4
    plot(t, X(:, i), '-', 'Color', cols(i, :), 'LineWidth', 1.5, 'DisplayName', sprintf('x_%d(t)', i - 1))
end
for i = 1 : 4
    plot(t, pi_(i) * ones(size(t)), ':', 'Color', cols(i, :), 'LineWidth', 1.5, 'HandleVisibility', 'off')
end
xlabel('t (hours)');
ylabel('probability');
title('Distribution Starting Empty')
legend();
box on;
hold off;

subplot(1, 2, 2);
semilogy(t, TV, 'k-', 'LineWidth', 1.5)
hold on;
semilogy(t, C * exp(lam2 * t), 'r--', 'LineWidth', 1.5)
semilogy(t, 0.01 * ones(size(t)), ':', 'LineWidth', 1.5)
semilogy([tstar tstar], [min(TV) 1], ':', 'LineWidth', 1.5)
plot(tstar, 0.01, 'ko', 'MarkerFaceColor', 'k')
xlabel('t (hours)')
ylabel('TV Distance (log scale)')
title('Total Variation Distance to \pi')
legend({'TV(t)', sprintf('C e^{\\lambda_2 t}, \\lambda_2 = %.3f', lam2)});
box on;
hold off;

print(fig, 'APPN4560P3', '-dpng', '-r300')