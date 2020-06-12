% This wrapper calls the NORST function to perform the task of Background
% Recovery for real data (videos) in presence of occlusions and corruptions

clear all;
clc;


addpath('YALL1_v1.4') % ell-1 minimization package
addpath('data') % contains of sample videos
addpath('PROPACK/');
addpath('export_fig-master/')
%% Loading video data matrix
video = ["Curtain","SwitchLight","Lobby"];
PATH = '../norst-miss/data';
load([PATH,'/',char(video(1)),'.mat'])

L = I;  % video with foreground
Train = DataTrain;  % video with background only
L = Train(:,101:end); 

%% Parameter Initialization
[n,m] = size(L);
r = 30;
t_max = m;
alpha = 60;

%% Generating missing entries' support
% (a)-Bernoulli Model:

rho = 0.1; %denotes fraction of missing entries
BernMat = rand(n, t_max);
T = 1 .* (BernMat <= 1 - rho);
    
% (b)-Moving Model: (rectangular object moving along the width of the video)

height = imSize(1);
width = imSize(2);    

% b-a = rectangle's width
a = floor(height/2) - 2;
b = a + 4;

% frames with moving object
% T = ones(n,m);
% 
% idx_frame = [width * 0 + 1 : width * floor(t_max/width)];
% 
% smin = 0;
% smax = 1;
% for j = idx_frame   
%     for i = smin:smax
%         T(height*i+ a : height*i + b ,j) = zeros(b-a+1,1);
%     end
%     smax = smax+1;
%     if(smax - smin > width/4)
%         smin = smin + 1;
%     end
% 
%     if(smax >= width)
%         smax = 1;
%         smin = 0;
%     end
% end

M = L .* T; % corrupted version of the video frames

%M = L;
%% Calling NORST
fprintf('NORST\n')

% algorithm parameters
K = 3;
ev_thresh = 2e-3;
omega = 15;
tol = 1e-3; % tolerance in cgls and ncrpca functions

% mean subtraction
mu = mean(Train(:,1:100),2);
M_norst = M - mu;

t_norst = tic;
% initialization of true subspace
fprintf('Initialization...\t');
%P_init = orth(ncrpca(L(:,1:100), r, tol, 100));
[P_init, ~] = svds(L(:,100), r);
%P_init = randn(n, r);
%fprintf('Subspace initialized\n');

[FG,BG] = NORST_video(M_norst, mu, T, P_init, ev_thresh, alpha, K, omega,tol);

t_NORST = toc(t_norst);                

%% Display the reconstructed video
%DisplayVideo(L, T, M, BG, imSize, 'recovered_BackGround_video_norst.avi')
