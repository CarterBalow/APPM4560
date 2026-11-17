clear; clc; close all;

lam = 12;
mu = 3;
rho = lam / mu;
T = 500;

t = 0;
n = 0;
max_states = 50;
state_time = zeros(1, max_states);

while t < T
    total_rate = lam + n * mu;
    dt = -log(rand()) / total_rate;

    if t + dt > T
        dt = T - t;
    end
    if n + 1 > length(state_time)
        state_time(n + 1) = 0;
    end

    state_time(n + 1) = state_time(n + 1) + dt;
    t = t + dt;

    if t < T
        if rand() < lam / total_rate
            n = n + 1;
        else
            n = n - 1;
        end
    end
end

last_idx = find(state_time > 0, 1, 'last');
states = 0 : (last_idx - 1);
emp_prob = state_time(1 : last_idx) / T;
theo_prob = (exp(-rho) .* (rho .^ states)) ./ factorial(states);

max_diff = max(abs(emp_prob - theo_prob));
time_avg = sum(states .* emp_prob);

c = 0;
while (1 - sum((exp(-rho) .* (rho .^ (0 : c))) ./ factorial(0 : c))) > 0.01
    c = c + 1;
end

fprintf('Largest prob diff: %.4f\n', max_diff);
fprintf('Time avg of X_t: %.4f\n', time_avg);
fprintf('Smallest fleet size c: %d (P(X > %d) = %.4f)\n', c, c, 1 - sum(1 - sum((exp(-rho) .* (rho .^ (0 : c))) ./ factorial(0 : c))));

fig = figure();
fig.Theme = 'light';
bar(states, [emp_prob', theo_prob']);
xlabel('Number of Bikes Out (n)')
ylabel('Probabilit / Time Fraction')
title('Gillespie Simulation vs. Poisson(4)')
legend('Simulated Fraction of Time', 'Theoretical Fraction of Time')
grid on;

print(fig, 'APPM4560HW10P2', '-dpng', '-r300')