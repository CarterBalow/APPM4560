clear; clc; close all;

N = 2e5;
T = 5;
q = 0.4;
r_star = 3;
r_sky = 2;
d_star = q * r_star;
d_sky = q * r_sky;

%% first way
star1 = zeros(N, 1);
sky1 = zeros(N, 1);

for win = 1 : N
    t = -log(rand) / r_star;
    while t <= T
        if rand < q
            star1(win) = star1(win) + 1;
        end
        t = t - log(rand) / r_star;
    end

    t = -log(rand) / r_sky;
    while t <= T
        if rand < q
            sky1(win) = sky1(win) + 1;
        end
        t = t - log(rand) / r_sky;
    end
end

%% second way
star2 = zeros(N, 1);
sky2 = zeros(N, 1);

for win = 1 : N
    t = -log(rand) / d_star;
    while t <= T
        star2(win) = star2(win) + 1;
        t = t - log(rand) / d_star;
    end

    t = -log(rand) / d_sky;
    while t <= T
        sky2(win) = sky2(win) + 1;
        t = t - log(rand) / d_sky;
    end
end

%% analysis
tot1 = star1 + sky1;
tot2 = star2 + sky2;

c1 = corrcoef(star1, sky1);
c2 = corrcoef(star2, sky2);

idx1 = (tot1 == 10);
idx2 = (tot2 == 10);

fprintf('Thin:   star %.4f / %.4f, sky %.4f / %.4f, corr = %.4f, total %.4f / %.4f, P(6|10) = %.4f\n', ...
    mean(star1), var(star1), mean(sky1), var(sky1), c1(1, 2), mean(tot1), var(tot1), mean(star1(idx1) == 6));
fprintf('Direct: star %.4f / %.4f, sky %.4f / %.4f, corr = %.4f, total %.4f / %.4f, P(6|10) = %.4f\n', ...
    mean(star2), var(star2), mean(sky2), var(sky2), c2(1, 2), mean(tot2), var(tot2), mean(star2(idx2) == 6));

%% detection times
Ts = -log(rand(N, 1)) / d_star;
Tk = -log(rand(N, 1)) / d_sky;

first = min(Ts, Tk);
is_star = (Ts < Tk);
M = max(Ts, Tk);

c_first = corrcoef(first, double(is_star));

fprintf('P(star first) = %.4f, P(star first, < 1 s) = %.4f, E[M] = %.4f, corr = %.4f\n', ...
    mean(is_star), mean(is_star & first < 1), mean(M), c_first(1, 2));