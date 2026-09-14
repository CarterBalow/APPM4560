clear; clc; close all;
%% 1
A = [0.9 0.1 0; 0 0.75 0.25; 0.5 0 0.5];
A_50 = A^50;

%% 2
p = [0.5 0.5 0 0 0 0; 1/3 0 1/3 0 1/3 0; 0 0 0.25 0.75 0 0; 0 0 1 0 0 0; 0 0 0 0 0 1; 0 0 0 0 1 0];
p_20 = p^20;
p_21 = p^21;

%% 4
% initializing
R = 2 * 10^4;
T = 1 * 10^4;
d0 = 10;

lamb_alive_b = ones(R, 1);
alive_count_b = zeros(T, 1);
lamb_alive_c = ones(R, 1);
alive_count_c = zeros(T, 1);

X1 = zeros(R, 1);
X2 = zeros(R, 1);
Y1 = ones(R, 1) * d0;
Y2_lion1 = ones(R, 1) * d0;
Y2_lion2 = ones(R, 1) * d0;

for t = 1 : T
    % part b (1 lion)
    dx1 = 2 * (rand(R, 1) > 0.5) - 1;
    X1 = X1 + dx1 .* lamb_alive_b;
    dy1 = 2 * (rand(R, 1) > 0.5) - 1;
    Y1 = Y1 + dy1 .* lamb_alive_b;

    lamb_alive_b = lamb_alive_b .* (X1 ~= Y1);
    alive_count_b(t) = sum(lamb_alive_b);

    % part c (2 lions)
    dx2 = 2 * (rand(R, 1) > 0.5) - 1;
    X2 = X2 + dx2 .* lamb_alive_c;
    dy2_lion1 = 2 * (rand(R, 1) > 0.5) - 1;
    Y2_lion1 = Y2_lion1 + dy2_lion1 .* lamb_alive_c;
    dy2_lion2 = 2 * (rand(R, 1) > 0.5) - 1;
    Y2_lion2 = Y2_lion2 + dy2_lion2 .* lamb_alive_c;

    lamb_alive_c = lamb_alive_c .* (X2 ~= Y2_lion1) .* (X2 ~= Y2_lion2);
    alive_count_c(t) = sum(lamb_alive_c);
end

S1 = alive_count_b / R;
S2 = alive_count_c / R;
S1_squared = S1 .^ 2;
t = (1 : T)';

fit = (t >= 10^2) & (t <= 10^4);
t_fit = t(fit);

poly_1 = polyfit(log(t_fit), log(S1(fit)), 1);
poly_2 = polyfit(log(t_fit), log(S2(fit)), 1);
beta_1 = -poly_1(1);
beta_2 = -poly_2(1);

S1_cont = erf(d0 ./ (2 * sqrt(t)));

% plotting
fig = figure();
fig.Theme = 'light';
hold on;
loglog(t, S1, 'b-', 'LineWidth', 1.5)
loglog(t, S1_cont, 'c--', 'LineWidth', 1.5)
loglog(t, S2, 'r-', 'LineWidth', 1.5);
loglog(t, S1_squared, 'm--', 'LineWidth', 1.5)

xlabel('Time step (t)')
ylabel('Survival Probability S_N(t)')
title(sprintf('Lamb Survival Rate (Fit Window: 10^2 \x2264 t \x2264 10^4)\n\\beta_1 = %.3f, \\beta_2 = %.3f', beta_1, beta_2))
legend('S_1(t) (Simulated)', 'S_1(t) (Continuum)', 'S_2(t) (Simulated)', 'S_1(t)^2', 'Location', 'best')
grid on;

print(fig, 'LambSurvivability', '-dpng', '-r300')