clear; clc; close all;

P = [0, 0.5, 0.5, 0, 0;
    0.25, 0, 0, 0.5, 0.25;
    0.75, 0, 0, 0, 0.25;
    0, 0, 0, 1, 0;
    0, 0, 0, 0, 1];
P_cum = cumsum(P, 2);

Q = P(1 : 3, 1 : 3);
R = P(1 : 3, 4 : 5);

state_names = {'U', 'I', 'M'};
exact_h = [0.5, 0.625, 0.375];
exact_g = [4, 2, 4];
exact_tauF = [4, 1.8, 5];
exact_tauA = [4, 7/3, 3.4];

num_sims = 10^4;
sim_data = struct();

for s_idx = 1 : 3
    endings = zeros(num_sims, 1);
    times = zeros(num_sims, 1);
    for i = 1 : num_sims
        curr_state = s_idx;
        steps = 0;
        while curr_state < 4
            steps = steps + 1;
            u = rand();
            curr_state = find(P_cum(curr_state, :) >= u, 1, 'first');
        end
        endings(i) = curr_state;
        times(i) = steps;
    end
    
    sim_data(s_idx).endings = endings;
    sim_data(s_idx).times = times;
end

tabulated = table();
for s_idx = 1:3
    endings = sim_data(s_idx).endings;
    times = sim_data(s_idx).times;
    
    h_sim = mean(endings == 4);
    g_sim = mean(times);
    
    tF = times(endings == 4);
    tA = times(endings == 5);
    
    tauF_sim = mean(tF);
    tauA_sim = mean(tA);
    
    tabulated = [tabulated; table(state_names(s_idx), exact_h(s_idx), h_sim, exact_g(s_idx), g_sim, exact_tauF(s_idx), tauF_sim, exact_tauA(s_idx), tauA_sim, ...
        'VariableNames', {'State', 'h_exact', 'h_sim', 'g_exact', 'g_sim', 'tauF_exact', 'tauF_sim', 'tauA_exact', 'tauA_sim'})];
end

disp(tabulated);

I_endings = sim_data(2).endings;
I_times = sim_data(2).times;

I_times_F = I_times(I_endings == 4);
I_times_A = I_times(I_endings == 5);

max_bins = 15;
n_vec = 1 : max_bins;

pmf_F = zeros(1, max_bins);
pmf_A = zeros(1, max_bins);

for n = 1 : max_bins
    if n == 1
        pmf_F(n) = R(2, 1) / exact_h(2);
        pmf_A(n) = R(2, 2) / (1 - exact_h(2));
    else
        Q_pow_R = (Q^(n - 1)) * R;
        pmf_F(n) = Q_pow_R(2, 1) / exact_h(2);
        pmf_A(n) = Q_pow_R(2, 2) / (1 - exact_h(2));
    end
end

fig = figure();
fig.Theme = 'light';

subplot(1, 2, 1);
hold on;
histogram(I_times_F, 'BinEdges', 0.5 : 1 : (max_bins + 0.5), 'Normalization', 'pdf');
stem(n_vec, pmf_F, 'r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
title('Conditioned on F (Start: I)');
xlabel('Steps to Absorption (n)');
ylabel('Probability Density');
legend('Simulated', 'Exact PMF');
xlim([0 max_bins + 1]);
grid on;

subplot(1, 2, 2);
hold on;
histogram(I_times_A, 'BinEdges', 0.5 : 1 : (max_bins + 0.5), 'Normalization', 'pdf');
stem(n_vec, pmf_A, 'r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
title('Conditioned on A (Start: I)');
xlabel('Steps to Absorption (n)');
ylabel('Probability Density');
legend('Simulated', 'Exact PMF');
xlim([0 max_bins + 1]);
grid on;

print(fig, 'ProteinFoldingChain', '-dpng', '-r300')