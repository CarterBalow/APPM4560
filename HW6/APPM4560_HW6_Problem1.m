clear; clc; close all;

for m = 3 : 12
    p = zeros(m, m);
    for i = 1 : m -1
        p(i, i + 1) = 1;
    end
    p(m, 1) = 0.5;
    p(m, 2) = 0.5;
end

B = p > 0;
M = 1;
while ~all(B(:))
    B = B * p > 0;
    M = M + 1;
end

lam = eig(p);
[sorted_abs, ~] = sort(abs(lam), 'descend');
slem = sorted_abs(2);

wielandt_bound = (m - 1) ^ 2 + 1;
