%% INITIALIZATION
clear; clc; close all;

N = 10^4;

lam1 = 0.5;
c1 = 4 / exp(1);
[samples1, rate1, mean_time1] = rejection_sampling(lam1, c1, N);

lam2 = 0.2;
c2 = 1 / (0.16 * exp(1));
[samples2, rate2, mean_time2] = rejection_sampling(lam2, c2, N);

%% PLOTTING
fig = figure();
fig.Theme = 'light';
x_vals = linspace(0, 10, 1000);
f_target = x_vals .* exp(-x_vals);

subplot(1, 2, 1)
hold on;
histogram(samples1, 50, 'Normalization', 'pdf');
plot(x_vals, f_target, 'r-', 'LineWidth', 1.5)
title('\lambda = 0.5')
xlabel('x');
ylabel('Density');
grid on;

subplot(1, 2, 2)
hold on;
histogram(samples2, 50, 'Normalization', 'pdf');
plot(x_vals, f_target, 'r-', 'LineWidth', 1.5)
title('\lambda = 0.2')
xlabel('x');
ylabel('Density');
grid on;

print(fig, 'RejectionSampling', '-dpng', '-r300')


%% REJECTION SAMPLING FUNC
function [samples, rate, mean_time] = rejection_sampling(lamda, c, N)
    samples = zeros(N, 1);
    count = 0;
    total_attempts = 0;

    tic;
    while count < N
        total_attempts = total_attempts + 1;

        U1 = rand();
        U2 = rand();
        X = -log(U1) / lamda;

        fX = X * exp(-X);
        gX = lamda * exp(-lamda * X);
        ratio = fX / (c * gX);

        if U2 < ratio
            count = count + 1;
            samples(count) = X;
        end
    end
    total_time = toc;
    
    rate = N / total_attempts;
    mean_time = total_time / N;
end