% % function script_Exp2
% Description:
%   main entrance for audio signal experiments
% ==============================================================

% Add path for codes
addpath(genpath(['.', filesep]));    

StrICA_lst = {'fastica','fasticagauss','jade','infomax','1fica', 'giica', 'pegi', 'pimd'};
Seed = 1:50;

StrICA_lst = {'fastica','fasticagauss','jade','1fica', 'giica', 'pegi', 'pimd'};
Seed = 1:5;

PMD = [0, 8, 10, 20, -8, -10, -20];
P_noise = [0.1,0.2:0.1:0.5];
PMD=sort(PMD);
P_noise=sort(P_noise);

% MAIN CALL
mainAudSig(PMD, StrICA_lst, P_noise, Seed);
