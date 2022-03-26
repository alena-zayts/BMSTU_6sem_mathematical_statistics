X = [7.76,6.34,5.11,7.62,8.84,4.68,8.65,6.90,8.79,6.61,6.62,7.13,6.75,7.28,7.74,7.08,5.57,8.20,7.78,7.92,6.00,4.88,6.75,6.56,7.48,8.51,9.06,6.94,6.93,7.79,5.71,5.93,6.81,5.76,5.88,7.05,7.22,6.67,5.59,6.57,7.28,6.22,6.31,5.51,6.69,7.12,7.40,6.86,7.28,6.82,7.08,7.52,6.81,7.55,4.89,5.48,7.74,5.10,8.17,7.67,7.07,5.80,6.10,7.15,7.88,9.06,6.85,4.88,6.74,8.76,8.53,6.72,7.21,7.42,8.29,8.56,9.25,6.63,7.49,6.67,6.79,5.19,8.20,7.97,8.64,7.36,6.72,5.90,5.53,6.44,7.35,5.18,8.25,5.68,6.29,6.69,6.08,7.42,7.10,7.14,7.10,6.60,6.35,5.99,6.17,9.05,6.01,7.77,6.27,5.81,7.80,9.89,4.39,6.83,6.53,8.15,6.68,6.87,6.31,6.83];

n = length(X);
fprintf('n = %i\n', n);


fprintf('\na)\n');
M_max = max(X);
M_min = min(X);
fprintf('M_max = %f\n', M_max);
fprintf('M_min = %f\n', M_min);


fprintf('\nб)\n');
R = M_max - M_min;
fprintf('R = %f\n', R);


fprintf('\nв)\n');
mu = mean(X);
s_sqr = var(X); % эквивалентно var(X, 0); =std(X)^2; Деленние на N-1
fprintf('mu = %f\n', mu);
fprintf('s_sqr = %f\n', s_sqr);


fprintf('\nг)\n');
m = floor(log2(n)) + 2;
delta = R / m;
fprintf('m = %i, delta = %f\n\n', m, delta);

ni_array = zeros([1 m]);
ai_array = zeros([1 m]);
bi_array = zeros([1 m]);
for i = 1:m-1
    ai_array(i) = M_min + (i - 1) * delta;
    bi_array(i) = M_min + i * delta;
    ni_array(i) = count_ni(i, ai_array(i), bi_array(i), X);
end
ai_array(m) = M_min + (m - 1) * delta;
bi_array(m) = M_max;
ni_array(m) = count_ni(m, ai_array(m), bi_array(i) + 1, X); %правая граница - включительно

fprintf('\n');
for i = 1:m-1
    fprintf('J%i = [%f; %f); n%i = %i\n', i, ai_array(i), bi_array(i), i, ni_array(i)); 
end
fprintf('J%i = [%f; %f]; n%i = %i\n', m, ai_array(m), bi_array(m), m, ni_array(m)); 


fprintf('\nд)\n');
figure;
hold on;
% гистограмма
h_array = ni_array * (1 / (n * delta));
ai_array(end+1) = M_max;
hist = histogram('BinEdges', ai_array,'BinCounts', h_array); 
%histogram(X, m) не производит деления на (n * delta); 
% график функции плотности распределения вероятностей нормальной случайной 
% величины с заданным математическим ожиданием и дисперсией 
sigma = sqrt(s_sqr);
x_array = (M_min - 1):(sigma / 100):(M_max + 1);
f = normpdf(x_array, mu, sigma);
plot(x_array, f, 'r', 'LineWidth', 2);
grid on;
xlabel('x');
ylabel('fn(x)');


fprintf('\ne)\n');
figure;
hold on;
% график эмпирической функции распределения
% можно стандартной функцией: ecdf_values = ecdf(X)
z_array = unique(X); % отсортированные уникальные значения из X
z_array(end + 1) = z_array(end) + 1; % чтобы график не обрывался
nt_array = zeros([1 length(z_array)]);
for z_index = 1:length(z_array)
    nt_array(z_index) = count_nt(z_array(z_index), X) / n;
end
stairs(z_array, nt_array, 'b');
% график функции распределения нормальной случайной 
% величины с заданным математическим ожиданием и дисперсией 
F = normcdf(x_array, mu, sigma); 
plot(x_array, F, 'r');
grid on;
xlabel('x');
ylabel('Fn(x)');



function [ni] = count_ni(i, a, b, X)
    tolerance = 1e-4;
    ni = 0;
    fprintf('J%i-ому интервалу принадлежат значения {', i);
    for x_index = 1:length(X)
        if (X(x_index) - a >= tolerance) && (X(x_index) - b < tolerance)
            fprintf('%f, ', X(x_index));
            ni = ni + 1;
        end
    end
    fprintf('}\n');
end

function [nt] = count_nt(t, X)
    tolerance = 1e-4;
    nt = 0;
    for x_index = 1:length(X)
        if (X(x_index) - t < tolerance)
            nt = nt + 1;
        end
    end
end
