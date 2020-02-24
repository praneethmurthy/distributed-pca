function P_hat = simpleEVD_fed(X, r, power_iter, taubatch, sig_c)
%   This MATLAB code implements the Simple-EVD algorithm and returns a
%   basis matrix for the new subspace
%
%   Input:
%   X = data matrix, and then compute empirical covariance
%   r = target rank of output
%
%   Output:
%   P_hat = basis matrix for output (m x r)

%[P_hat, ~] = svds(X, r);

n = size(X, 1);
u_noise = randn(n, r);
%conv_noise = zeros(1, power_iter+1);
%conv_noise(1) = sin(subspace(u_true, u_noise));
for ii = 1 : power_iter
    ch_noise = sig_c * randn(n, r);
    u_noise = X * (X' * u_noise) + ch_noise;
    
    %conv_noise(ii+1) = sin(subspace(u_true, u_noise));
    if(~mod(ii, taubatch))
        [u_noise,~] = qr(u_noise, 0);
    end
end
P_hat = u_noise;
end
