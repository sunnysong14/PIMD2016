function [S_est, W_est, A_est, C] = mainICA_fun(strICA, pMD, X, Noise_cov, A)
% Usage: 
%   All ICA methods. Centering and whitening are conducted within each of
%   them for this noisyICA work. 
% 
% Inputs
%   strICA -- String name of ICA method.
%   X -- original mixing data. It will become $X_pn$, the p-level gaussian
%   noise corrupted X in noisy case.
% 
%   %--------The followings are unique for opt-SINR---------------------%
%   Noise_cov -- covariance matrix of noise. Used, e.g. for opt-SINR.
%   A -- original mixing matrix
% Outputs
%   S_est -- Estimated sources.
%   W_est -- Inverse of mixing matrix A, demixing matrix. W_est*X=S.

%%
% Inner product matrix C for PEGI and PIMD:
C = -1;     % invalid
A_est = -1; % invalid

% [core-1] Mean Displace
X_pMD = X - pMD*(mean(X,2)*ones(1,size(X,2)));

% The global variable TAG_CENTER is set to zero here, in order to control 
% and remove centering step for all investigated ICA methods.
global TAG_CENTER
TAG_CENTER = 0;

% [core-2] ICA methods
switch lower(strICA)
    case 'fastica' % 'pow3' (default)
        [S_est,~,W_est] = fastica_LY(X_pMD);
    case 'fasticagauss' % 'gauss': the robust version [Hyvarinen, TNN'99].
        [S_est,~,W_est] = fastica_LY(X_pMD,'g','gauss');
    case 'fasticatanh' % 'tanh'
        [S_est,~,W_est] = fastica_LY(X_pMD,'g','tanh');
    case 'jade'
        W_est = jadeR_LY(X_pMD);
        S_est = W_est*X_pMD;
    case 'infomax'
        W_est = infomaxICA_LY(X_pMD);
        S_est = W_est*X_pMD;
    case '1fica'
        W_est = FICA1_LY(X_pMD);
        S_est = W_est*X_pMD;
    case 'giica' 
        [S_est,W_est,~,A_est,C] = GIICA_LY(X_pMD,'quasi-orthogonalize','SINR variant',0,'verbose',0);
    case 'pegi'
        [S_est,W_est,~,A_est,C] = GIICA_LY(X_pMD,'pseudo-Euclidean IP','SINR variant',1,'verbose',0);
    case 'pimd' %Our method
        [S_est, W_est, ~, A_est, C] = GIICA_LY(X_pMD, 'real-Euclidean IP', 'SINR variant', 1, 'verbose', 0);
    case 'ainv'
        W_est = eye(size(A)) / A;
        S_est = W_est*X;
    case 'sinropt' % Voss-nips15 ICAExperimentRuns().
        W_est = A' / (Noise_cov+A*A');
        S_est = W_est*X;        
    otherwise
        error('Error: No such "caseICA" defined.');
end%END OF SWITCH
end%END OF FUNCTION
