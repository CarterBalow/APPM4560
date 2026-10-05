clear; clc; close all;

%% (b)
N = 4;
states = 0 : N;

% a = 1 = b
P1 = zeros(N + 1, N + 1);
for k = 0 : N
    idx = k + 1;
    if k < N
        P1(idx, idx + 1) = (4 - k) / 4;
    end
    if k > 0
        P1(idx, idx - 1) = k / 4;
    end
    P1(idx, idx) = 1 - sum(P1(idx, :));
end

n_max = 60;
q = zeros(n_max + 1, N + 1);
q0 = [1; 0; 0; 0; 0];
q(1, :) = q0;

for n = 1 : n_max
    q(n + 1, :) = q(n, :) * P1;
end

q50 = q(51, :);
q51 = q(52, :);

qn_2 = q(:, 3);
running_avqn2 = zeros(n_max + 1, 1);
for n = 0 : n_max
    running_avqn2(n + 1) = mean(qn_2(1: n + 1));
end

fig1 = figure();
fig1.Theme = 'light';
n_vals = 0 : n_max;
hold on;
plot(n_vals, qn_2, 'LineWidth', 1.5)
plot(n_vals, q(:, 5), 'LineWidth', 1.5)
plot(n_vals, running_avqn2, 'LineWidth', 1.5)
title('Potassium Channel Periodic Dynamics (a = 1 = b)');
xlabel('Step (n)');
ylabel('Probability');
legend('q_n(2)', 'q_n(4)', 'Running Average q_n(2)', 'Location', 'best');
grid on;

%% (c)

a = 0.3;
b = 0.1;

P2 = zeros(N + 1, N + 1);
for k = 0 : N
    idx = k + 1;
    if k < N
        P2(idx, idx + 1) = a * (4 - k) / 4;
    end
    if k > 0
        P2(idx, idx - 1) = b * k / 4;
    end
    P2(idx, idx) = 1 - sum(P2(idx, :));
end

A = P2' - eye(N + 1);
A(end + 1, :) = ones(1, N + 1);
temp = zeros(N + 2, 1);
temp(end) = 1;
pi_final = (A \ temp)';

n_max_qq = 300;
qq = zeros(n_max_qq + 1, N + 1);
qq(1, :) = q0;

n_smallest = -1;
for n = 1 : n_max_qq
    qq(n + 1, :) = qq(n, :) * P2;
    max_diff = max(abs(qq(n + 1, :) - pi_final));
    if max_diff < 1E-6 && n_smallest == -1
        n_smallest = n;
    end
end

sorted_eigs = sort(abs(eig(P2)), 'descend');
slem = sorted_eigs(2);

fig2 = figure();
fig2.Theme = 'light';
n_vals2 = 0 : n_max_qq;
hold on;
plot(n_vals2, qq, 'LineWidth', 1.5)
title('Potassium Channel Aperiodic Dynamics (a = 0.3, b = 0.1)');
xlabel('Step (n)');
ylabel('Probability');
legend('q_n(0)', 'q_n(1)', 'q_n(2)', 'q_n(3)', 'q_n(4)', 'Location', 'best');
grid on;

print(fig1, 'PeriodicChannel', '-dpng', '-r300')
print(fig2, 'AperiodicChannel', '-dpng', '-r300')