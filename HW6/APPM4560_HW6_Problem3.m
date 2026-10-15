clear; clc; close all;

%% (c)
R0_vals = linspace(0.5, 3, 500);
h1_100 = zeros(size(R0_vals));
h1_1000 = zeros(size(R0_vals));

N1 = 100;
K1 = 10;
for i = 1 : length(R0_vals)
    R0 = R0_vals(i);
    rho = zeros(1, K1);
    rho(1) = 1;

    for j = 1 : K1 - 1
        rho(j + 1) = rho(j) * (N1 / (R0 * (N1 - j)));
    end
    h1_100(i) = 1 / sum(rho);
end

N2 = 1000;
K2 = 100;
for i = 1 : length(R0_vals)
    R0 = R0_vals(i);
    rho = zeros(1, K2);
    rho(1) = 1;

    for j = 1 : K2 - 1
        rho(j + 1) = rho(j) * (N2 / (R0 * (N2 - j)));
    end
    h1_1000(i) = 1 / sum(rho);
end

asymptote = max(0, 1 - 1 ./ R0_vals);

fig1 = figure();
fig1.Theme = 'light';
hold on;
grid on;
plot(R0_vals, h1_100)
plot(R0_vals, h1_1000)
plot(R0_vals, asymptote)
xlabel('R_0')
ylabel('h_1')
title('Breakout Probability')
legend('N = 100, K = 10', 'N = 1000, K = 100', '1 - 1/R_0', 'Location', 'best')

print(fig1, 'BreakoutProb(c)', '-dpng', '-r300')

%% (d)

gam = 1;
dt = 0.005;
R0s = [0.5, 0.8, 1.5, 3.0];
e1dt = zeros(size(R0s));
num = 1;

for x = R0s
    beta = x * gam / N1;
    size_sys = K1 - 1;

    A = zeros(size_sys, size_sys);
    B = ones(size_sys, 1);
    for n = 1 : size_sys
        pn = beta * n * (N1 - n) * dt;
        qn = gam * n * dt;
        A(n, n) = pn + qn;
        
        if n > 1
            A(n, n - 1) = -qn;
        end
        if n < size_sys
            A(n, n + 1) = -pn;
        end
    end

    e = A \ B;
    e1dt(num) = e(1) * dt;
    num = num + 1;
end

R0_sim = 0.5;
beta_sim = R0_sim * gam / N1;
num_sims = 100000;
T_steps = zeros(num_sims, 1);

for i = 1 : num_sims
    X = 1;
    steps = 0;
    while X > 0 && X < K1
        pn = beta_sim * X * (N1 - X) * dt;
        qn = gam * X * dt;
        r = rand();
        if r < pn
            X = X + 1;
        elseif r < pn + qn
            X = X - 1;
        end
        steps = steps + 1;
    end
    T_steps(i) = steps;
end

sim_e1dt = mean(T_steps) * dt;

fprintf('Act e1dt = %.5f, Sim e1dt = %.5f', e1dt(1), sim_e1dt)