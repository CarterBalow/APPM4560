clear; clc; close all;

data = readtable('big12.csv');

all_teams = unique([data.home; data.away]);
m = length(all_teams);

team_names = all_teams;

games2025 = data(data.season == 2025, :);
games2026 = data(data.season == 2026, :);


%% 2025 data
pi2025 = 0.5 * ones(m, m);
h_sum25 = zeros(m, m);
counts25 = zeros(m, m);

for game = 1 : height(games2025)
    id_h = find(strcmp(team_names, games2025.home{game}));
    id_a = find(strcmp(team_names, games2025.away{game}));
    S = games2025.home_pts(game);
    T = games2025.away_pts(game);

    x_h = (S + 1) / (S + T + 2);
    hij_h = 0.5 + 0.5 * sign(x_h - 0.5) * sqrt(abs(2 * x_h - 1));
    x_a = (T + 1) / (S + T + 2);
    hij_a = 0.5 + 0.5 * sign(x_a - 0.5) * sqrt(abs(2 * x_a - 1));

    h_sum25(id_h, id_a) = h_sum25(id_h, id_a) + hij_h;
    h_sum25(id_a, id_h) = h_sum25(id_a, id_h) + hij_a;
    counts25(id_h, id_a) = counts25(id_h, id_a) + 1;
    counts25(id_a, id_h) = counts25(id_a, id_h) + 1;
end

for i = 1 : m
    for j = 1 : m
        if counts25(i, j) > 0
            pi2025(i, j) = h_sum25(i, j) / counts25(i, j);
        end
    end
end


%% 2026 data
w = 1;
played2026 = games2026(games2026.completed == 1, :);

n = zeros(m, 1);
h_sum26 = zeros(m, m);

for game = 1 : height(played2026)
    id_h = find(strcmp(team_names, played2026.home{game}));
    id_a = find(strcmp(team_names, played2026.away{game}));
    S = played2026.home_pts(game);
    T = played2026.away_pts(game);

    x_h = (S + 1) / (S + T + 2);
    hij_h = 0.5 + 0.5 * sign(x_h - 0.5) * sqrt(abs(2 * x_h - 1));
    x_a = (T + 1) / (S + T + 2);
    hij_a = 0.5 + 0.5 * sign(x_a - 0.5) * sqrt(abs(2 * x_a - 1));

    h_sum26(id_h, id_a) = h_sum26(id_h, id_a) + hij_h;
    h_sum26(id_a, id_h) = h_sum26(id_a, id_h) + hij_a;
    n(id_h) = n(id_h) + 1;
    n(id_a) = n(id_a) + 1;
end


%% Prior combos
pi_flat = 0.5 * ones(m, m);

si = [15; 12; 8; 14; 10; 9; 11; 7; 13; 6; 5; 4; 16; 3; 2; 1];
pi_poll = zeros(m, m);
for i = 1 : m
    for j = 1 : m
        pi_poll(i, j) = si(i) / (si(i) + si(j));
    end
end

A_flat = zeros(m, m);
A_2025 = zeros(m, m);
A_poll = zeros(m, m);

for i = 1 : m
    denom = n(i) + w * (m - 1);
    for j = 1 : m
        if i == j
            continue;
        end
        A_flat(i, j) = (h_sum26(i, j) + w * pi_flat(i, j)) / denom;
        A_2025(i, j) = (h_sum26(i, j) + w * pi2025(i, j)) / denom;
        A_poll(i, j) = (h_sum26(i, j) + w * pi_poll(i, j)) / denom;
    end
end


%% Analysis
[V1, D1] = eig(A_flat);
eigenvalues1 = diag(D1);
[~, idx1] = sort(abs(eigenvalues1), 'descend');
l_flat = real(eigenvalues1(idx1(1)));
ratio_flat = abs(eigenvalues1(idx1(2)) / l_flat);
r_flat = real(V1(:, idx1(1)));
r_flat = r_flat / sum(r_flat);

[V2, D2] = eig(A_2025);
eigenvalues2 = diag(D2);
[~, idx2] = sort(abs(eigenvalues2), 'descend');
l_2025 = real(eigenvalues2(idx2(1)));
ratio_2025 = abs(eigenvalues2(idx2(2)) / l_2025);
r_2025 = real(V2(:, idx2(1)));
r_2025 = r_2025 / sum(r_2025);

[V3, D3] = eig(A_poll);
eigenvalues3 = diag(D3);
[~, idx3] = sort(abs(eigenvalues3), 'descend');
l_poll = real(eigenvalues3(idx3(1)));
ratio_poll = abs(eigenvalues3(idx3(2)) / l_poll);
r_poll = real(V3(:, idx3(1)));
r_poll = r_poll / sum(r_poll);


%% Plotting
fprintf('Prior 1 (Flat): lambda* = %.4f, |lambda2/lambda*| = %.4f\n', l_flat, ratio_flat);
fprintf('Prior 2 (2025): lambda* = %.4f, |lambda2/lambda*| = %.4f\n', l_2025, ratio_2025);
fprintf('Prior 3 (Poll): lambda* = %.4f, |lambda2/lambda*| = %.4f\n', l_poll, ratio_poll);

fig = figure();
fig.Theme = 'light';
bar_data = [r_flat, r_2025, r_poll];
bar(bar_data);
set(gca, 'XTick', 1:m, 'XTickLabel', team_names, 'XTickLabelRotation', 45);
ylabel('r');
title('Keener Ratings');
legend({'Flat', '2025', 'Poll'}, 'Location', 'Best');
grid on;


%% Match predictions
target_games = games2026(games2026.completed == 0, :);
if ~isempty(target_games)
    for game = 1 : height(target_games)
        home_t = target_games.home{game};
        away_t = target_games.away{game};
        id_h = find(strcmp(team_names, home_t));
        id_a = find(strcmp(team_names, away_t));
        
        if r_flat(id_h) >= r_flat(id_a)
            fprintf('Week %d: %s defeats %s\n', target_games.week(game), home_t, away_t);
        else
            fprintf('Week %d: %s defeats %s\n', target_games.week(game), away_t, home_t);
        end
    end
end
