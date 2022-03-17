X = [7.76,6.34,5.11,7.62,8.84,4.68,8.65,6.90,8.79,6.61,6.62,7.13,6.75,7.28,7.74,7.08,5.57,8.20,7.78,7.92,6.00,4.88,6.75,6.56,7.48,8.51,9.06,6.94,6.93,7.79,5.71,5.93,6.81,5.76,5.88,7.05,7.22,6.67,5.59,6.57,7.28,6.22,6.31,5.51,6.69,7.12,7.40,6.86,7.28,6.82,7.08,7.52,6.81,7.55,4.89,5.48,7.74,5.10,8.17,7.67,7.07,5.80,6.10,7.15,7.88,9.06,6.85,4.88,6.74,8.76,8.53,6.72,7.21,7.42,8.29,8.56,9.25,6.63,7.49,6.67,6.79,5.19,8.20,7.97,8.64,7.36,6.72,5.90,5.53,6.44,7.35,5.18,8.25,5.68,6.29,6.69,6.08,7.42,7.10,7.14,7.10,6.60,6.35,5.99,6.17,9.05,6.01,7.77,6.27,5.81,7.80,9.89,4.39,6.83,6.53,8.15,6.68,6.87,6.31,6.83];

% а
M_max = max(X);
M_min = min(X);

% б
R = M_max - M_min;

% в
MX = mean(X);
DX = var(X); % sigma == std == sqrt(var(arg))

% г
m = floor(log2(length(X))) + 2;
h = histogram(X, m);
%disp(h);

% д
sigma = std(X);
x = (M_min - 1):(sigma / 100):(M_max + 1);
f = normpdf(x, MX, sigma); % normal probability distribution function
figure;
heights = h.Values / (sum(h.Values) * h.BinWidth);
centers = [];
for i = 1:(length(h.BinEdges) - 1)
    centers = [centers, (h.BinEdges(i + 1) + h.BinEdges(i)) / 2];
end
%disp(centers);
hold on;
bar(centers, heights, 1); % ширина относительная:)
plot(x, f, 'g', 'LineWidth', 2);

% е)
F = normcdf(x, MX, sigma); % normal cumulative distribution function
figure;
hold on;
ecdf(X); % empiric cumulative distribution function
plot(x, F, 'r');
