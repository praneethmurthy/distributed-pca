%% some codes to test out the channel noise sort of setting for power method

%%move all the codes that try to generate more noise along the true
%%direction to a different file...



clear
clc
close all
rng shuffle

%% define data model parameters
n = 1000;
r = 50;
signal_energy = 1.1;
noise_energy = 1.0;
ch_noise_energy = 1e-8;


%% generate the "rectangular" data points, and also the sample covariance
u_orth = orth(randn(n, 2 * (r+1)));
u_true = u_orth(:,1:r);
Y = u_orth(:,1:r+1) * diag([repmat(signal_energy, 1, r), noise_energy]) * u_orth(:, r+2:end)';
X = Y * Y';

%X =  X / (0.5 * signal_energy^2);
%verifying that the sample PC is close to the true PC
[u_init, s_init, v_init] = svds(X, r);
fprintf('SE after adding small noise: %d \n', ...
    sin(subspace(u_true, u_init)))

figure;
plot(svd(X))
title('singular values of original matrix')


%% vanilla power method -- with normalization
power_iter = 500;
taubatch = 10;
u_vanilla_norm = randn(n, r);
u_vanilla_norm = u_vanilla_norm / norm(u_vanilla_norm);
conv_vanilla_norm = zeros(1, power_iter+1);
conv_vanilla_norm(1) = sin(subspace(u_init, u_vanilla_norm));
for ii = 1 : power_iter
    u_vanilla_norm = X * u_vanilla_norm;
    [u_vanilla_norm, ~ ] = qr(u_vanilla_norm, 0);
    %u_vanilla_norm = u_vanilla_norm/norm(u_vanilla_norm);
    conv_vanilla_norm(ii+1) = sin(subspace(u_true, u_vanilla_norm));
end 

figure;
subplot(221)
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
    if(~mod(ii, taubatch))
        %[u_vanilla,~] = qr(u_vanilla, 0);
        u_vanilla = orth(u_vanilla);
    end
    conv_vanilla(ii+1) = sin(subspace(u_true, u_vanilla));
end 

subplot(222)
plot([1: power_iter+1], log10(conv_vanilla))
axis tight
title('Convergence vanilla PM -- without norm')
stry = '$$\log(SE(\hat{u}_t, u))$$';
strx = '$$\mathrm{power\ iterations} (t) $$';
ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18) 
xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18) 

fprintf('SE for vanilla power method (without norm): %d \n', ...
    conv_vanilla_norm(end))


%% channel noise power method -- with normalization
u_noise_norm = randn(n, r);
%u_noise_norm = u_noise_norm / norm(u_noise_norm);
conv_noise_norm = zeros(1, power_iter+1);
conv_noise_norm(1) = sin(subspace(u_true, u_noise_norm));
for ii = 1 : power_iter
    ch_noise = ch_noise_energy * randn(n, r);
    u_noise_norm = X * u_noise_norm + ch_noise;
    %u_noise_norm = u_noise_norm/norm(u_noise_norm);
    [u_noise_norm,~] = qr(u_noise_norm, 0);
    conv_noise_norm(ii+1) = sin(subspace(u_true, u_noise_norm));
end 

subplot(223)
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
    u_noise = X * u_noise + ch_noise;
    
    conv_noise(ii+1) = sin(subspace(u_true, u_noise));
    if(~mod(ii, taubatch))
        norm(u_noise)
        %[u_noise,~] = qr(u_noise, 0);
        u_noise = orth(u_noise);
    end
end 

subplot(224)
plot([1: power_iter+1], log10(conv_noise))
axis tight
title('with channel noise but no norm')
stry = '$$\log(SE(\hat{u}_t, u))$$';
strx = '$$\mathrm{power\ iterations} (t) $$';
ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18) 
xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18) 

fprintf('SE for noisy power method (without norm): %d \n', conv_noise(end))



% figure;
% subplot(211)
% 
% plot([1: power_iter+1], log10(conv_vanilla_norm), 'kd', 'LineStyle', '--', 'MarkerSize', 6, 'LineWidth', 2)
% hold
% plot([1: power_iter+1], log10(conv_noise_norm),'rs', 'LineStyle', '-.', 'MarkerSize', 6, 'LineWidth', 2)
% %plot([1: power_iter+1], log10(conv_noise_sig_norm), 'go', 'LineStyle', '-', 'MarkerSize', 6, 'LineWidth', 2)
% axis tight
% grid on
% % legend('noiseless', 'with channel noise', 'with biased channel noise')
% l1 = legend('noiseless', 'with channel noise');
% l1.FontSize = 15;
% stry = '$$\log(SE({u}_t, u_1^*))$$';
% strx = '$$\mathrm{power\ iterations} (t) $$';
% ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18) 
% xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18) 
% title('PM with normalization')
% 
% subplot(212)
% plot([1: power_iter+1], log10(conv_vanilla), 'k+', 'LineStyle', '--', 'MarkerSize', 6, 'LineWidth', 2)
% hold
% plot([1: power_iter+1], log10(conv_noise),'rs', 'LineStyle', '-.', 'MarkerSize', 6, 'LineWidth', 2)
% %plot([1: power_iter+1], log10(conv_noise_sig), 'go', 'LineStyle', '-', 'MarkerSize', 6, 'LineWidth', 2)
% axis tight
% grid on
% % legend('noiseless', 'with channel noise', 'with biased channel noise')
% l1 = legend('noiseless', 'with channel noise');
% l1.FontSize = 15;
% 
% stry = '$$\log(SE({u}_t, u_1^*))$$';
% strx = '$$\mathrm{power\ iterations} (t) $$';
% ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18) 
% xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18) 
% title('PM without normalization')


figure;
plot(linspace(1, power_iter+1, (power_iter)/taubatch+1), log10(conv_noise_norm(1:taubatch:end)),'k+', 'LineStyle', '--', 'MarkerSize', 6, 'LineWidth', 2)
hold
plot(linspace(1, power_iter+1, (power_iter)/taubatch+1), log10(conv_noise(1:taubatch:end)),'rs', 'LineStyle', '-.', 'MarkerSize', 6, 'LineWidth', 2)
%plot([1: power_iter+1], log10(conv_noise_sig_norm), 'go', 'LineStyle', '-', 'MarkerSize', 6, 'LineWidth', 2)
axis tight
grid on
% legend('noiseless', 'with channel noise', 'with biased channel noise')
l1 = legend('normalize every iter', 'normalize after taubatch');
l1.FontSize = 15;
stry = '$$\log(SE({u}_t, u_1^*))$$';
strx = '$$\mathrm{power\ iterations} (t) $$';
ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18) 
xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18) 
title('PM with normalization')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the algorithm works even with uniform r.v. initialization as opposed to
%gaussian init. This is interesting, and should ideally be used to
%strengtehn the theory



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%some old codes regarding biased power method!
% val = 1e-3;
% tmpmatrix = val * (u_true * u_true') + sqrt((1 - val^2))/n * (eye(n) - u_true * u_true');
% trace(tmpmatrix)
% R = chol(tmpmatrix);


% %% ``biased channel noise'' with norm
% u_noise_sig_norm = randn(n, r);
% u_noise_sig_norm = u_noise_sig_norm / norm(u_noise_sig_norm);
% conv_noise_sig_norm = zeros(1, power_iter+1);
% conv_noise_sig_norm(1) = sin(subspace(u_true, u_noise_sig_norm));
% for ii = 1 : power_iter
%     %ch_noise_sig_norm = ch_noise_energy * randn(n, r);
%     ch_noise_sig_norm = ch_noise_energy * R * randn(n, r);
%     ch_noise_sig_norm = ch_noise_sig_norm + val * u_true;
%     %abs(ch_noise_sig_norm' * u_true)
%     %ch_noise_sig_norm = (val * (u_true * u_true') + ...
%     %    sqrt(1 - val^2) * (eye(n) - u_true * u_true')) * ch_noise_sig_norm;
%     u_noise_sig_norm = X * u_noise_sig_norm + ch_noise_sig_norm;
%     u_noise_sig_norm = u_noise_sig_norm/norm(u_noise_sig_norm);
%     conv_noise_sig_norm(ii+1) = sin(subspace(u_true, u_noise_sig_norm));
% end 
% 
% subplot(325)
% plot([1: power_iter+1], log10(conv_noise_sig_norm))
% axis tight
% title('biased channel noise with norm')
% stry = '$$\log(SE(\hat{u}_t, u))$$';
% strx = '$$\mathrm{power\ iterations} (t) $$';
% ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18) 
% xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18) 
% 
% fprintf('SE for biased noisy power method (without norm): %d \n', conv_noise_sig_norm(end))
% 
% 
% 
% %% ``biased channel noise''
% u_noise_sig = randn(n, r);
% conv_noise_sig = zeros(1, power_iter+1);
% conv_noise_sig(1) = sin(subspace(u_true, u_noise_sig));
% for ii = 1 : power_iter
%     ch_noise_sig = ch_noise_energy * R * randn(n, r);
%     ch_noise_sig = ch_noise_sig + val * u_true;
%     %ch_noise_sig = (val * (u_true * u_true') + ...
%     %     sqrt(1 - val^2) * (eye(n) - u_true * u_true')) * ch_noise_sig + ch_noise_sig;
%     %abs(ch_noise_sig_norm' * u_true)
%     u_noise_sig = X * u_noise_sig + ch_noise_sig;
%     conv_noise_sig(ii+1) = sin(subspace(u_true, u_noise_sig));
% end 
% 
% subplot(326)
% plot([1: power_iter+1], log10(conv_noise_sig))
% axis tight
% title('biased channel noise but no norm')
% stry = '$$\log(SE(\hat{u}_t, u))$$';
% strx = '$$\mathrm{power\ iterations} (t) $$';
% ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18) 
% xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18) 
% 
% fprintf('SE for biased noisy power method (without norm): %d \n', conv_noise_sig(end))
%how to fix the rank r case
