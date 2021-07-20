clear;
clc;

addpath('NORST-rmc-master/');
addpath('NORST-rmc-master/PROPACK');



%% Algorithms to run
NORST = 1;
NORST_OFFLINE = 0;

%% Parameter Initialization
n = 1000;
t_max = 3000;
alpha = 60;
f = 100;
MC = 1;
t_calc_pca = [1:alpha: t_max];
t_calc = t_calc_pca;

%NORST
temp_SE_NORST = zeros(length(t_calc_pca), MC);
temp_err_L_NORST = zeros(t_max, MC);

temp_SE_NORST_fed = zeros(length(t_calc_pca), MC);
temp_err_L_NORST_fed = zeros(t_max, MC);


temp_SE_NORST_base = zeros(length(t_calc_pca), MC);
temp_err_L_NORST_base = zeros(t_max, MC);

t_NORST = 0;
err_L_fro_NORST = zeros(MC,1);

for mc = 1 : MC
     fprintf('Monte-Carlo iteration %d in progress \n', mc);

    %% Generating support set using: a. Moving Object Model    b. Bernoulli Model
    % a. Moving Object Model
%     s = 100;            
%     S = zeros(n, t_max);
%     rho_s = 1;
%     b0 = 0.1;
%     beta = ceil(b0 * alpha1);
%     x_max = 25;
%     x_min = 15;
%     alpha1 = 100;
%     num_changes = floor((t_max -t_train)/beta);
% 
%     num_changes1 = min(floor(alpha1 / beta), ceil(n/s));
% 
%     flag = 0;
%     ii1 = 1;
%     fval1 = 0;
%     for ii = 1 : num_changes
%         if(~flag)   %%downward motion
%             if(ii1 <= num_changes1)
%                 bind = fval1 + (ii1 - 1) * s/rho_s + 1;
%                 sind = min(bind - 1 + s, n);
%                 ii1 = ii1 + 1;
%                 if(ii1 == num_changes1 + 1)
%                     flag = 1;
%                     ii1 = 1;
%                     fval2 = bind;
%                 end
%             end
%         else
%             if(ii1 <= num_changes1)
%                 bind = max(fval2 - (ii1 - 1) * s/rho_s , 1);
%                 sind = bind - 1 + s;
%                 ii1 = ii1 + 1;
%                 if(ii1 == num_changes1 + 1)
%                     flag = 0;
%                     ii1 = 1;
%                 end
%             end
%         end
%         idx = bind : sind;
%         jdx = (ii-1) * beta + 1 : ii * beta;
%         S(idx, jdx) = x_min + ...
%             (x_max - x_min) * rand(length(idx), beta);
%         T(idx, jdx) = 1;
%     end
%     
    % b. Bernoulli Model
    rho = 0.1; % fraction of missing entries
    BernMat = rand(n, t_max);
    T = 1 .* (BernMat <= 1 - rho); % observed entries' support 
 
    %% Generating low-rank matrix
    r_0 = 10;
    L = zeros(n, t_max);
    
    lambda_min = sqrt(f)/2;
    lambda_max = sqrt(f);
    
    offset = 0; %if offset is not zero, eigenvalues (Lambda) are varying in time;
    diag_entries1 = offset +  [linspace(lambda_max, lambda_min, r_0)];
    diag_entries2 = -offset + [linspace(lambda_max, lambda_min, r_0)];

    coeff_train = zeros(r_0, t_max);
    for cc = 1 : r_0
        coeff_train(cc, 1:2:end-1) = -diag_entries1(cc) + ...
            2 * diag_entries1(cc) * rand(1, t_max/2);
        
        coeff_train(cc, 2:2:end) = -diag_entries2(cc) + ...
            2 * diag_entries2(cc) * rand(1, t_max/2);
    end
    
    P = orth(randn(n, r_0));
    delta_t = 100;
    U0 = P;
    subspace_size = 60;
    U_track = cell(ceil(t_max/subspace_size),1);
    for i=1:length(U_track)
        Btemp1 = randn(n);
        B1 = (Btemp1 - Btemp1')/2;
        
        t_1 = (i-1)*subspace_size + 1;
        t_2 = min(i*subspace_size,t_max);
        
        U_track{i} = U0;
        L(:, t_1:t_2) = U0 * coeff_train(:,t_1:t_2);
        
        U = expm(delta_t*B1)*U0;        
        U0 = U;
    end 
    
%     P = orth(randn(n, r_0));
%     delta_t = 1e-4;
%     U0 = P;
%     U_track = cell(t_max,1);
%     Btemp1 = randn(n);
%     B1 = (Btemp1 - Btemp1')/2;
%     
%     rot_mat = expm(delta_t*B1);
%     
%     tmp_se = zeros(size(U_track));
%     for i=1:length(U_track)
%         
%         
%         
%         U_track{i} = U0;
%         L(:, i) = U0 * coeff_train(:,i);
%         U = rot_mat * U0;        
%         %tmp_se(i) = Calc_SubspaceError(U, U0);
%         U0 = U;
%     end 
    
    
    eps_noise = 0; % noise
    L = L + eps_noise * (rand(n,t_max) - 0.5);
    M = L .* T ;
    
    %% Algorithm parameters for NORST
    if(NORST == 1)
    fprintf('\tNORST\t');
    K = 25;
    ev_thresh = 7.5961e-04;
    tol = 1e-16;
    overlap_step = alpha; % if it is set to alpha then windows don't overlap
    R = 0; % number of reuse 
    
    P_init = zeros(n,r_0);
    
   % (M,mu,T_obs,P_init,ev_thresh,alpha,K,omega,tol)
    
    t_norst_fed = tic;
    [L_hat_fed, S_hat_fed, t_hat_fed, P_track_full_fed] =  ...
        NORST_nodet(M, T, P_init, ev_thresh, alpha, K, 1, tol);
    t_norst_fed = toc(t_norst_fed);
    
    t_norst = tic;
    [L_hat, P_hat, S_hat, t_hat, P_track_full, t_calc] =  ...
        NORST_random(M, T, r_0, ev_thresh, alpha, K,R,overlap_step);
    t_NORST = toc(t_norst)
    err_L_fro_NORST(mc) = norm(L-L_hat,'fro')/norm(L,'fro');
    
    t_norst_base = tic;
    [L_hat_base, S_hat_base, t_hat_base, P_track_full_base] =  ...
        baseline_PCA(M, T, P_init, ev_thresh, alpha, K, 1, tol);
    t_norst_base = toc(t_norst_base);
    
    end

    %% Compute Performance Metrics
    %frobenius norm errors
if(NORST == 1)
    temp_err_L_NORST(:, mc) = sqrt(mean((L - L_hat).^2, 1)) ...
        ./ sqrt(mean(L.^2, 1));
    temp_err_L_NORST_fed(:, mc) = sqrt(mean((L - L_hat_fed).^2, 1)) ...
        ./ sqrt(mean(L.^2, 1));
end

    %subspace errors
   for jj = 1 : length(t_calc_pca)
            %tt = t_calc_pca(jj);
            tt = ceil(t_calc_pca(jj)/subspace_size);
            if(NORST == 1)           
            temp_SE_NORST(jj, mc) = ...
                Calc_SubspaceError(P_track_full{jj}, U_track{tt});
            temp_SE_NORST_fed(jj, mc) = ...
                Calc_SubspaceError(P_track_full_fed{jj}, U_track{tt});
            temp_SE_NORST_base(jj, mc) = ...
                Calc_SubspaceError(P_track_full_base{jj}, U_track{tt});
            end
            
            if(NORST_OFFLINE == 1)           
            temp_SE_NORST_off(jj, mc) = ...
                Calc_SubspaceError(P_track_full{jj}, U_track{tt});
            end
    end
fprintf('\n')
end

err_SE_NORST = mean(temp_SE_NORST, 2);
err_L_NORST = mean(temp_err_L_NORST, 2);

err_SE_NORST_fed = mean(temp_SE_NORST_fed, 2);
err_L_NORST_fed = mean(temp_err_L_NORST_fed, 2);

err_SE_NORST_base = mean(temp_SE_NORST_base, 2);
err_L_NORST_base = mean(temp_err_L_NORST_base, 2);

figure
strx = 'j';
stry = '$$\log_{10} dist(\hat{P}_{j}, P_{j})$$';

semilogy(t_calc_pca(1:1:end),err_SE_NORST(1:1:end),'-*r','LineWidth',2,'MarkerSize',10);
hold
semilogy(t_calc_pca(1:1:end),err_SE_NORST_fed(1:1:end),'-*b','LineWidth',2,'MarkerSize',10);
semilogy(t_calc_pca(1:1:end),err_SE_NORST_base(1:1:end),'-*k','LineWidth',2,'MarkerSize',10);

grid on

axis tight
xlabel(strx, 'Interpreter', 'LaTeX', 'FontSize', 20);
ylabel(stry, 'Interpreter', 'LaTeX', 'FontSize', 20);
legend('NORST', 'NORST-nodet', 'Baseline-PCA')