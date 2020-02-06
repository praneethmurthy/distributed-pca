%% some codes to test out the channel noise sort of setting for power method

%%move all the codes that try to generate more noise along the true
%%direction to a different file...


%%dependence on taubatch


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
power_iter = 500;
%taubatch = 10;

taubatchrange = [5, 10, 20, 50];
%taubatchrange = [2, 4, 6, 8];%5, 10, 20, 50];
%taubatchrange = [5, 10];

MC = 50;
conv_noise_mat = zeros(power_iter +1, length(taubatchrange), MC);



for mc = 1 : MC
    fprintf('%d MC  iteration\n', mc);
    %% generate the "rectangular" data points, and also the sample covariance
    u_orth = orth(randn(n, 2 * (r+1)));
    u_true = u_orth(:,1:r);
     Y = u_orth(:,1:r+1) * diag([repmat(signal_energy, 1, r), noise_energy]) * u_orth(:, r+2:end)';
    %Y = u_orth(:,1:r+1) * diag([linspace(signal_energy, 10 * signal_energy, r), noise_energy]) * u_orth(:, r+2:end)';
    X = Y * Y';
    
    %X =  X / (0.5 * signal_energy^2);
    %verifying that the sample PC is close to the true PC
    [u_init, s_init, v_init] = svds(X, r);
    fprintf('SE after adding small noise: %d \n', ...
        sin(subspace(u_true, u_init)))
    
%     figure;
%     plot(svd(X))
%     title('singular values of original matrix')
    
    
    %% channel noise power method -- normalization every iteration
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
    
    fprintf('SE for noisy power method: %d \n', conv_noise_norm(end))
    
    %% channel noise power method -- normalization only for certain iterations
    
    ctr = 1;
    for taubatch = taubatchrange
        u_noise = randn(n, r);
        conv_noise = zeros(1, power_iter+1);
        conv_noise(1) = sin(subspace(u_true, u_noise));
        for ii = 1 : power_iter
            ch_noise = ch_noise_energy * randn(n, r);
            u_noise = X * u_noise + ch_noise;
            
            conv_noise(ii+1) = sin(subspace(u_true, u_noise));
            if(~mod(ii, taubatch))
                [u_noise,~] = qr(u_noise, 0);
            end
        end
        conv_noise_mat(:, ctr, mc) = conv_noise';
        ctr = ctr + 1;
    end
end


fprintf('SE for noisy power method (without norm): %d \n', conv_noise(end))


conv_noise_mat = mean(conv_noise_mat, 3);


figure;
plot(linspace(1, power_iter+1, (power_iter)/taubatchrange(1)+1), log10(conv_noise_norm(1:taubatchrange(1):end)),'k+', 'LineStyle', '--', 'MarkerSize', 6, 'LineWidth', 2)
hold
plot(linspace(1, power_iter+1, (power_iter)/taubatchrange(1)+1), log10(conv_noise_mat(1:taubatchrange(1):end,1)),'rs', 'LineStyle', '-.', 'MarkerSize', 6, 'LineWidth', 2)
plot(linspace(1, power_iter+1, (power_iter)/taubatchrange(2)+1), log10(conv_noise_mat(1:taubatchrange(2):end,2)),'b+', 'LineStyle', '-.', 'MarkerSize', 6, 'LineWidth', 2)
plot(linspace(1, power_iter+1, (power_iter)/taubatchrange(3)+1), log10(conv_noise_mat(1:taubatchrange(3):end,3)),'gd', 'LineStyle', '-.', 'MarkerSize', 6, 'LineWidth', 2)
plot(linspace(1, power_iter+1, (power_iter)/taubatchrange(4)+1), log10(conv_noise_mat(1:taubatchrange(4):end,4)),'mo', 'LineStyle', '-.', 'MarkerSize', 6, 'LineWidth', 2)
%plot([1: power_iter+1], log10(conv_noise_sig_norm), 'go', 'LineStyle', '-', 'MarkerSize', 6, 'LineWidth', 2)
axis tight
grid on
% legend('noiseless', 'with channel noise', 'with biased channel noise')
l1 = legend('normalize every iter', 'taubatch = 5', 'taubatch = 10', 'taubatch=20', 'taubatch=50');
l1.FontSize = 15;
stry = '$$\log(SE({\hat U}_t, U))$$';
strx = '$$\mathrm{power\ iterations} (t) $$';
ylabel(stry, 'Interpreter', 'latex', 'FontSize', 18)
xlabel(strx, 'Interpreter', 'latex', 'FontSize', 18)
title('PM with normalization')
