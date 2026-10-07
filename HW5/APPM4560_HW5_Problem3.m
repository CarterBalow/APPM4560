clear; clc; close all;

%% (a)
p = [0 0.5 0.5 0 0 0 0 0;
     0 0 0.5 0.5 0 0 0 0;
     0.5 0 0 0 0.5 0 0 0;
     1/3 0 1/3 0 1/3 0 0 0;
     0 1/3 0 0 0 1/3 1/3 0;
     0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 1 0];

q0 = (1/8) * ones(1, 8);

q100 = q0 * (p^100);
q101 = q0 * (p^101);

%% (b)

d = 0.85;
G = d * p + ((1 - d) / 8) * ones(8, 8);

q_curr = q0;
iter = 0;
while true
    q_next = q_curr * G;
    iter = iter + 1;
    if norm(q_next - q_curr, 1) < 1E-10
        pi_final = q_next;
        break;
    end
    q_curr = q_next;
end

A = G' - eye(8);
A(end + 1, :) = ones(1, 8);
temp = [zeros(8, 1); 1];
pi_solve = (A \ temp)';

%% (c)

R = 100;
T = 10 ^ 5;
surfers = zeros(R, T);
surfers(:, 1) = 1;

% simulation 1 (surfer 1 tracking)
G_cum = cumsum(G, 2);

for t = 1 : T - 1
    curr_pages = surfers(:, t);
    uniforms = rand(R, 1);
    next_pages = sum(uniforms > G_cum(curr_pages, :), 2) + 1;
    surfers(:, t + 1) = next_pages;
end

surfer1_counts = histcounts(surfers(1, :), 1 : 9);
pi_hat_s1 = surfer1_counts / T;

% simulation 2 (error analysis)
T_vals = round(logspace(2, 5, 100));
rms_errs = zeros(size(T_vals));

for i = 1 : length(T_vals)
    T_curr = T_vals(i);
    all_ests = zeros(R, 8);

    for r = 1 : R
        counts = histcounts(surfers(r, 1 : T_curr), 1 : 9);
        all_ests(r, :) = counts / T_curr;
    end

    max_errs = max(abs(all_ests - pi_final), [], 2);
    rms_errs(i) = sqrt(mean(max_errs .^ 2));
end

logT = log10(T_vals');
logErr = log10(rms_errs');
coeffs = polyfit(logT, logErr, 1);

% plotting
fig1 = figure();
fig1.Theme = 'light';
bar([pi_final; pi_hat_s1]')
xlabel('Page Value')
ylabel('Probability')
title('True $\pi$ vs. Surfer 1 Estimate $\hat{\pi}$', 'Interpreter', 'latex')
legend('True \pi', 'Surfer 1 Estimate')
grid on;

fig2 = figure();
fig2.Theme = 'light';
hold on;
loglog(T_vals, rms_errs, 'b-', 'LineWidth', 1.5);
loglog(T_vals, 10 ^ coeffs(2) * T_vals .^ coeffs(1), 'r--', 'LineWidth', 1.5);
xlabel('Time Steps (T)')
ylabel('RMS Error')
title(['Error Decay vs. Time Steps (Fitted Slope: ', num2str(coeffs(1)), ')'])
legend('Simulated RMS Error', 'Fitted Slope')
grid on;
ax = gca; 
ax.XScale = 'log';
ax.YScale = 'log';

print(fig1, 'PageSurferSimulation', '-dpng', '-r300')
print(fig2, 'PageSurfingRMSError', '-dpng', '-r300')