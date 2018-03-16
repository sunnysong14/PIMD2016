function [X_pn, Noise_cov] = injctGauNoise(X, p_noise, seed)
% Usage:
%   Inject additive Gaussian noise into X.
% 	The covariance of noise is assumed to be p_noise*I where p_noise=sigma^2.
% 
% INPUT/OUTPUT PARAMETERS:
%   X -- Observer matrix. Each row: a source, and each column: a sampling.
%   p_noise -- a positive constant defining the noise powers with its value between 0
% and 1. In math, p_noise=sigma^2.
%   X_pn: Corrupted X with p_noise gaussian noise.
%%
rng(seed); % control random

[nSrc, nSmp] = size(X);
Noise_mean = zeros(nSrc, nSmp);

% Cov(noise)
Noise_cov = p_noise*eye(nSrc);

% Inject Gaussian noise
X_pn = X + mvnrnd(Noise_mean', Noise_cov)';
end  % END OF FUNCTION
