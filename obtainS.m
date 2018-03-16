function S = obtainS(seed)
% Usage: 
%   Load audio signals.
%% #seed < 100
seed = rem(seed, 101);

% load source signals
flnm_source = '31audio14S.mat';
load(flnm_source); %>> 'S' & 'Info_audio14S'

% permute S's rows
rng(seed); %control random
Ind = randperm(size(S,1));
S = S(Ind, :);
%music show: j=3; sig=S(j,:); Fs=16000; player=audioplayer(sig,Fs); play(player)

% For now, only study nS=4
nS = 4;
S = S(1:nS,:);
end%END OF FUNCTION
