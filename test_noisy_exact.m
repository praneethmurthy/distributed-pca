%% some codes to test out the channel noise sort of setting for power method


clear
clc
%close all
rng shuffle

n = 100;
r = 1;

signal_energy = 5;
noise_energy = 1;
ch_noise_energy = 1;

u_orth = orth(randn(n, 2 * (r+1)));

u_true = u_orth(:,1);

val = 1e-3;


tmpmatrix = val * (u_true * u_true') + sqrt((1 - val^2))/n * (eye(n) - u_true * u_true');
trace(tmpmatrix)
R = chol(tmpmatrix);

%noise_temp = noise_energy * sqrt(1/n) * randn(n);

Y = u_orth(:,1:2) * diag([signal_energy, noise_energy]) * u_orth(:, 3:4)';

X = Y * Y';

%X = signal_energy^2 * (u_true * u_true') + (noise_temp + noise_temp');
%X = X/norm(X);

[u_init, s_init, v_init] = svds(X, r);

fprintf('SE after adding small noise: %d \n', ...
    sin(subspace(u_true, u_init)))

figure;
plot(svd(X))
title('singular values of original matrix')


%% vanilla power method -- with normalization
power_iter = 100;
u_vanilla_norm = randn(n, r);
u_vanilla_norm = u_vanilla_norm / norm(u_vanilla_norm);
conv_vanilla_norm = zeros(1, power_iter+1);
conv_vanilla_norm(1) = sin(subspace(u_init, u_vanilla_norm));
for ii = 1 : power_iter
    u_vanilla_norm = X * u_vanilla_norm;
    u_vanilla_norm = u_vanilla_norm/norm(u_vanilla_norm);
    conv_vanilla_norm(ii+1) = sin(subspace(u_true, u_vanilla_norm));
end 

figure;
subplot(321)
plot([1: power_iter+1], log10(conv_vanilla_norm))
axis tight
title('Convergence for vanilla PM')
stry = '$$\log(SE(\hat{u}_t, u))$$';
strx = '$$\mathrm{power\ iterations} (t) $$';
ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18) 
xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18) 

fprintf('SE for vanilla power method: %d \n', conv_vanilla_norm(end))

%% vanilla power method -- without normalization
u_vanilla = randn(n, r);
conv_vanilla = zeros(1, power_iter+1);
conv_vanilla(1) = sin(subspace(u_true, u_vanilla));
for ii = 1 : power_iter
    u_vanilla = X * u_vanilla;
    conv_vanilla(ii+1) = sin(subspace(u_true, u_vanilla));
end 

subplot(322)
plot([1: power_iter+1], log10(conv_vanilla))
axis tight
title('Convergence vanilla PM -- without norm')
stry = '$$\log(SE(\hat{u}_t, u))$$';
strx = '$$\mathrm{power\ iterations} (t) $$';
ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18) 
xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18) 

fprintf('SE for vanilla power method (without norm): %d \n', conv_vanilla_norm(end))


%% channel noise power method -- with normalization
u_noise_norm = randn(n, r);
u_noise_norm = u_noise_norm / norm(u_noise_norm);
conv_noise_norm = zeros(1, power_iter+1);
conv_noise_norm(1) = sin(subspace(u_true, u_noise_norm));
for ii = 1 : power_iter
    ch_noise = ch_noise_energy * randn(n, r);
    u_noise_norm = X * u_noise_norm + ch_noise;
    u_noise_norm = u_noise_norm/norm(u_noise_norm);
    conv_noise_norm(ii+1) = sin(subspace(u_true, u_noise_norm));
end 

subplot(323)
plot([1: power_iter+1], log10(conv_noise_norm))
axis tight
title('with channel noise and normalization')
stry = '$$\log(SE(\hat{u}_t, u))$$';
strx = '$$\mathrm{power\ iterations} (t) $$';
ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18) 
xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18) 

fprintf('SE for noisy power method: %d \n', conv_noise_norm(end))


%% channel noise power method -- without normalization

u_noise = randn(n, r);
conv_noise = zeros(1, power_iter+1);
conv_noise(1) = sin(subspace(u_true, u_noise));
for ii = 1 : power_iter
    ch_noise = ch_noise_energy * randn(n, r);
    %abs(ch_noise' * u_true)
    u_noise = X * u_noise + ch_noise;
    conv_noise(ii+1) = sin(subspace(u_true, u_noise));
end 

subplot(324)
plot([1: power_iter+1], log10(conv_noise))
axis tight
title('with channel noise but no norm')
stry = '$$\log(SE(\hat{u}_t, u))$$';
strx = '$$\mathrm{power\ iterations} (t) $$';
ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18) 
xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18) 

fprintf('SE for noisy power method (without norm): %d \n', conv_noise(end))

%% ``biased channel noise'' with norm
u_noise_sig_norm = randn(n, r);
u_noise_sig_norm = u_noise_sig_norm / norm(u_noise_sig_norm);
conv_noise_sig_norm = zeros(1, power_iter+1);
conv_noise_sig_norm(1) = sin(subspace(u_true, u_noise_sig_norm));
for ii = 1 : power_iter
    %ch_noise_sig_norm = ch_noise_energy * randn(n, r);
    ch_noise_sig_norm = ch_noise_energy * R * randn(n, r);
    ch_noise_sig_norm = ch_noise_sig_norm + val * u_true;
    %abs(ch_noise_sig_norm' * u_true)
    %ch_noise_sig_norm = (val * (u_true * u_true') + ...
    %    sqrt(1 - val^2) * (eye(n) - u_true * u_true')) * ch_noise_sig_norm;
    u_noise_sig_norm = X * u_noise_sig_norm + ch_noise_sig_norm;
    u_noise_sig_norm = u_noise_sig_norm/norm(u_noise_sig_norm);
    conv_noise_sig_norm(ii+1) = sin(subspace(u_true, u_noise_sig_norm));
end 

subplot(325)
plot([1: power_iter+1], log10(conv_noise_sig_norm))
axis tight
title('biased channel noise with norm')
stry = '$$\log(SE(\hat{u}_t, u))$$';
strx = '$$\mathrm{power\ iterations} (t) $$';
ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18) 
xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18) 

fprintf('SE for biased noisy power method (without norm): %d \n', conv_noise_sig_norm(end))



%% ``biased channel noise''
u_noise_sig = randn(n, r);
conv_noise_sig = zeros(1, power_iter+1);
conv_noise_sig(1) = sin(subspace(u_true, u_noise_sig));
for ii = 1 : power_iter
    ch_noise_sig = ch_noise_energy * R * randn(n, r);
    ch_noise_sig = ch_noise_sig + val * u_true;
    %ch_noise_sig = (val * (u_true * u_true') + ...
    %     sqrt(1 - val^2) * (eye(n) - u_true * u_true')) * ch_noise_sig + ch_noise_sig;
    %abs(ch_noise_sig_norm' * u_true)
    u_noise_sig = X * u_noise_sig + ch_noise_sig;
    conv_noise_sig(ii+1) = sin(subspace(u_true, u_noise_sig));
end 

subplot(326)
plot([1: power_iter+1], log10(conv_noise_sig))
axis tight
title('biased channel noise but no norm')
stry = '$$\log(SE(\hat{u}_t, u))$$';
strx = '$$\mathrm{power\ iterations} (t) $$';
ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18) 
xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18) 

fprintf('SE for biased noisy power method (without norm): %d \n', conv_noise_sig(end))



figure;
subplot(211)

plot([1: power_iter+1], log10(conv_vanilla_norm), 'kd', 'LineStyle', '--', 'MarkerSize', 6, 'LineWidth', 2)
hold
plot([1: power_iter+1], log10(conv_noise_norm),'rs', 'LineStyle', '-.', 'MarkerSize', 6, 'LineWidth', 2)
plot([1: power_iter+1], log10(conv_noise_sig_norm), 'go', 'LineStyle', '-', 'MarkerSize', 6, 'LineWidth', 2)
axis tight
grid on
legend('noiseless', 'with channel noise', 'with biased channel noise')
stry = '$$\log(SE(\hat{u}_t, u))$$';
strx = '$$\mathrm{power\ iterations} (t) $$';
ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18) 
xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18) 
title('with normalization')

subplot(212)
plot([1: power_iter+1], log10(conv_vanilla), 'k+', 'LineStyle', '--', 'MarkerSize', 6, 'LineWidth', 2)
hold
plot([1: power_iter+1], log10(conv_noise),'rs', 'LineStyle', '-.', 'MarkerSize', 6, 'LineWidth', 2)
plot([1: power_iter+1], log10(conv_noise_sig), 'go', 'LineStyle', '-', 'MarkerSize', 6, 'LineWidth', 2)
axis tight
grid on
legend('noiseless', 'with channel noise', 'with biased channel noise')
stry = '$$\log(SE(\hat{u}_t, u))$$';
strx = '$$\mathrm{power\ iterations} (t) $$';
ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18) 
xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18) 
title('without normalization')
