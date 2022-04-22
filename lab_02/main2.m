X = [7.76,6.34,5.11,7.62,8.84,4.68,8.65,6.90,8.79,6.61,6.62,7.13,6.75,7.28,...
7.74,7.08,5.57,8.20,7.78,7.92,6.00,4.88,6.75,6.56,7.48,8.51,9.06,6.94,6.93,...
7.79,5.71,5.93,6.81,5.76,5.88,7.05,7.22,6.67,5.59,6.57,7.28,6.22,6.31,5.51,...
6.69,7.12,7.40,6.86,7.28,6.82,7.08,7.52,6.81,7.55,4.89,5.48,7.74,5.10,8.17,...
7.67,7.07,5.80,6.10,7.15,7.88,9.06,6.85,4.88,6.74,8.76,8.53,6.72,7.21,7.42,...
8.29,8.56,9.25,6.63,7.49,6.67,6.79,5.19,8.20,7.97,8.64,7.36,6.72,5.90,5.53,...
6.44,7.35,5.18,8.25,5.68,6.29,6.69,6.08,7.42,7.10,7.14,7.10,6.60,6.35,5.99,...
6.17,9.05,6.01,7.77,6.27,5.81,7.80,9.89,4.39,6.83,6.53,8.15,6.68,6.87,6.31,...
6.83];

% X=[-2.54,-0.79,-4.27,-3.09,-3.82,-0.61,-0.64,-1.24,-1.73,-2.91,-1.48,-1.28,-0.37,-1.88,-2.19,-1.61,-1.52,-3.17,-1.36,-3.08,-3.11,-3.07,-1.57,-1.51,-2.37,-0.58,-3.05,-2.93,-1.01,-1.40,-2.06,-3.05,-1.84,-1.24,-1.89,-2.06,-1.59,-2.83,-1.07,-2.96,-3.17,-3.08,-0.49,-3.11,-3.14,-2.30,-3.99,-1.56,-1.28,-3.46,-2.63,-0.82,-2.18,-0.89,-3.08,-1.13,-1.62,-1.06,-2.98,-1.55,-1.49,-1.65,-1.45,-2.29,-0.85,-1.44,-2.87,-2.40,-2.13,-3.52,-1.42,-3.64,-3.47,-2.05,-2.39,-2.07,-0.80,-1.52,-3.92,-2.22,-0.78,-2.60,-1.78,-1.61,-1.65,-2.06,-3.33,-3.41,-1.97,-1.74,-2.04,0.01,-1.37,-3.15,-2.35,-3.66,-1.79,-2.56,-1.87,-1.06,-0.64,-2.49,-1.85,-1.40,-0.86,-0.17,-0.62,-2.85,-2.12,-1.17,-2.48,-1.65,-3.74,-2.87,-3.15,-1.89,-1.34,-4.33,-0.96,-1.79];


gamma = 0.9;
N = length(X);
% prompt = "Input gamma: ";
% gamma = input(prompt);

% 2
mu = find_mu(X);
s_sqr = find_s_sqr(X); 
mu_low = find_mu_low(mu, s_sqr, N, gamma);
mu_high = find_mu_high(mu, s_sqr, N, gamma);
s_sqr_low = find_sigma_sqr_low(s_sqr, N, gamma);
s_sqr_high = find_sigma_sqr_high(s_sqr, N, gamma);

fprintf("For given sample (N = % i):\n", N);
fprintf('mu = %f; s_sqr = %f\n', mu, s_sqr);
fprintf('mu_low = %f; mu_high = %f\n', mu_low, mu_high);
fprintf('s_sqr_low = %f; s_sqr_high = %f\n', s_sqr_low, s_sqr_high);


% 3
n_array = zeros([1 N]);
mu_array = zeros([1 N]);
mu_low_array = zeros([1 N]);
mu_high_array = zeros([1 N]);
s_sqr_array = zeros([1 N]);
s_sqr_low_array = zeros([1 N]);
s_sqr_high_array = zeros([1 N]);
for i = 1:N
    n_array(i) = i;

    mu_i = find_mu(X(1:i));
    mu_array(i) = mu_i;
    s_sqr_i = find_s_sqr(X(1:i));
    s_sqr_array(i) = s_sqr_i;

    mu_low_array(i) = find_mu_low(mu_i, s_sqr_i, i, gamma);
    mu_high_array(i) = find_mu_high(mu_i, s_sqr_i, i, gamma);
    
    s_sqr_low_array(i) = find_sigma_sqr_low(s_sqr_i, i, gamma);
    s_sqr_high_array(i) = find_sigma_sqr_high(s_sqr_i, i, gamma);
end

mu_const = mu * ones(N);
s_sqr_const = s_sqr * ones(N);

% a
plot(n_array, mu_const, n_array, mu_array, ...
    n_array, mu_low_array, n_array, mu_high_array);
xlabel('n');
ylabel('y');
xlim([1 N]);
legend('$\hat \mu(\vec x_N)$', '$\hat \mu(\vec x_n)$', ...
    '$\underline{\mu}(\vec x_n)$', '$\overline{\mu}(\vec x_n)$', ...
    'Interpreter', 'latex', 'FontSize', 14);
figure;

% b
plot(n_array, s_sqr_const, n_array, s_sqr_array, ...
    n_array, s_sqr_low_array, n_array, s_sqr_high_array);
xlabel('n');
ylabel('z');
xlim([1 N]);
legend('$\hat S^2(\vec x_N)$', '$\hat S^2(\vec x_n)$', ...
   '$\underline{\sigma}^2(\vec x_n)$', '$\overline{\sigma}^2(\vec x_n)$', ...
   'Interpreter', 'latex', 'FontSize', 14);

% functions
function [mu] = find_mu(X)
    mu = mean(X);
end

function [s_sqr] = find_s_sqr(X)
    s_sqr = var(X);
end

% tinv(a, n) - квантиль уровня a распределения Стьюдента с n степенями свободы.
function [mu_low] = find_mu_low(mu, s_sqr, n, gamma)
    mu_low = mu - sqrt(s_sqr) * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
end

function [mu_high] = find_mu_high(mu, s_sqr, n, gamma)
    mu_high = mu + sqrt(s_sqr) * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
end

%chi2inv(a, n) - квантиль уровня a распределения хи квадрат с n степенями свободы.
function [sigma_sqr_low] = find_sigma_sqr_low(s_sqr, n, gamma)
    sigma_sqr_low = (n - 1) * s_sqr / chi2inv((1 + gamma) / 2, n - 1);
end

function [sigma_sqr_high] = find_sigma_sqr_high(s_sqr, n, gamma)
    sigma_sqr_high = (n - 1) * s_sqr / chi2inv((1 - gamma) / 2, n - 1);
end
