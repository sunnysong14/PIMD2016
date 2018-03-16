function [X_pn, X, S, A, Noise_cov] = obtainICA_signal(p_noise, seed)
% Usage:
%   Obtain real-world ICA data: sources, mixing matrix, noise-free X and
%   noisy X.
% 
% Parameters:
%   nSmp -- #sampling of each source.
%   p_noise -- a positive constant defining the noise powers with its value between 0
% and 1, in formula it is p_noise=\sigma^2.
%   seed -- the seed with which we replicate the randomness.
% 
%   X_pn -- corrupted X with Gaussian noise level $p$.
%   X -- noise-free observation
%   S -- sources
%   A -- original mixing matrix
%   Noise_cov -- covariance matrix of noise. Used, e.g. for opt-SINR.

%% 
% load source
S = obtainS(seed);
nS = size(S,1);

% normalization
for ii = 1:size(S,1)
    Si = S(ii,:);
    S(ii,:) = Si./std(Si);
end

% mixing matrix
A = genA(nS,seed);

% noise-free X
X = A*S;

% inject noise to X
if p_noise > 0 
     [X_pn,Noise_cov] = injctGauNoise(X,p_noise,seed);
elseif p_noise == 0 
    X_pn = X;
    Noise_cov = 0;
else
    error(['Error in pNoise:',num2str(p_noise),'.\n']);
end
end % END OF FUNCTION

%% Test
%{
clear,clc, close all

caseScen = 5;
p_noise = 0.1; 
seed = 1; 

[X_pn, X, S, A, Noise_cov] = obtainICA_IS(caseScen,p_noise, seed);
%}
