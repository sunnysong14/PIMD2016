% [08/11/2016] for PIMD'acml16
% Liyan updates fastICA with $pMD$. 
% To change the codes as less as possible, I use Global variable.
% 
function [W2, W1]=FICA1_LY(x, ini, g)
% LY records:
%   This is a variant of the well known FastICA that is proposed for BSS in
%   offline (block processing) setup and a noisy environment. 
%   This algorithm combines Symmetric FastICA with test of saddle points to 
%   achieve fast global convergence and one unit refinement to obtain high 
%   noise rejection ability. It's computational complexity is similar to 
%   that of the original FastICA, i.e. it is fast.
%   url: http://itakura.ite.tul.cz/zbynek/1fica.htm
% 
% Liyan for NIPS16 04-14-2016, 05-12-2016: update centering and
% whitening, for PIMD-acml'16 08-09-2016: fix the randomness of 1FICA.
% 
% ========================================================================
% copyright: Zbynek Koldovsky, Petr Tichavsky
% The Best Blind MMSE Estimator
% 
% Abstract: This is a variant of the well known FastICA algorithm that is proposed to be used for blind source separation in off-line (block processing) setup and a noisy environment. The algorithm combines Symmetric FastICA with test of saddle points to achieve fast global convergence and one-unit refinement to obtain high noise rejection ability. It's computational complexity is similar to that of the original FastICA, i.e. it is fast. It is shown that the bias of the algorithm due to additive noise is reduced; it is proportional to ?^3, where ?^2 is the variance of the additive noise, while the bias of the other methods (namely of all methods using the orthogonality constraint, and even of recently proposed EFICA algorithm) is asymptotically proportional to ?^2. Thanks to the reduced bias, the novel algorithm exhibits a significantly lower symbol-error rate when it is applied to blindly separate mixtures of finite alphabet signals that are typical for communication systems.
% Corresponding papers:
%     Z. Koldovský and P. Tichavský, "Blind Instantaneous Noisy Mixture Separation with Best Interference-plus-noise Rejection", Proceedings of 7th International Conference on Independent Component Analysis (ICA2007), pp. 730-737, Sept. 2007. (here)
%     Z. Koldovský and P. Tichavský, "Asymptotic Analysis of Bias of FastICA-based Algorithms in Presence of Additive Noise", technical report nr. 2181, ÚTIA, AV ?R, Jan 2007.
% 
% Matlab codes: 1FICA for real signals; for complex signals FicaCPLX available at personal pages of Petr Tichavský
% Copyright: Zbyn?k Koldovský, Petr Tichavský
% version: 1.0  release: 19.1.2007

epsilon=0.000001;
[dim N]=size(x);
MaxIt=100;

repeat=1;
rot2d=[1/sqrt(2) 1/sqrt(2);-1/sqrt(2) 1/sqrt(2)];
SaddleTest=true;

if nargin<3
    g='tanh';
end

%<--[LY, 08-09-2016 for PIMD-acml'16]: randomness happens here.
% set random
iseed = 1;
rng(iseed);
% 
if nargin<2
    ini=randn(dim,dim);
end

test_of_saddle_points_nonln='tanh';

%% [08/11/2016] LY for Mean Displacement
% Note: TAG_CENTER is a global variable indicating centering or
% non-centering. The global variable is set outside in (e.g.) mainICA().

% Centering
global TAG_CENTER
if TAG_CENTER == 1
    for j=1:dim
        x(j,:)=x(j,:)-mean(x(j,:));
    end 
end

%% preprocessing, LY: whitening
C = cov(x');

%whitened data
x = C^(-1/2)*x;

%%% FastICA Symmetric approach with the test of saddle points
W=ini;
W=W*real(inv(W'*W)^(1/2));
while repeat
    crit=0; NumIt=0;
    while (1-min(crit)>epsilon && NumIt<MaxIt)% && sum(double(changed>10))<2)
        Wold=W;
        switch g
            case {'tanh','biga'}
                hypTan = tanh(x'*W);
                W=x*hypTan/N-ones(dim,1)*sum(1-hypTan.^2).*W/N;
            case 'pow3'
                W=(x*(pwr(x'*W,3)))/N-3*W;
            case 'gaus'
                U=x'*W;
                Usquared=U.^2;
                ex=exp(-Usquared/2);
                gauss=U.*ex;
                dGauss=(1-Usquared).*ex;
                W=x*gauss/N-ones(dim,1)*sum(dGauss).*W/N;
        end
        W=W*real(inv(W'*W)^(1/2));
        crit=abs(sum(W.*Wold));
        NumIt=NumIt+1;
    end %while iteration
    repeat=0;
    %%%The test of saddle points
    if SaddleTest
        SaddleTest=false; %%The test could be done only one times
        u=x'*W;
        switch test_of_saddle_points_nonln
            case 'tanh'
                table1=(mean(log(cosh(u)))-0.37456).^2;
            case 'gaus'
                table1=(mean(ex)-1/sqrt(2)).^2;
            case 'pow3'
                table1=(mean((pwr(u,4)))-3).^2;
        end
        rotated=zeros(1,dim);
        checked=1:dim; 
        for i=checked
            for j=checked(checked>i)
                if (~rotated(i) && ~rotated(j))
                    h=[u(:,i) u(:,j)]*rot2d;
                    switch test_of_saddle_points_nonln
                        case 'tanh'
                            ctrl=(mean(log(cosh(h)))-0.37456).^2;
                        case 'gaus'
                            ctrl=(mean(exp(-h.^2/2)-1/sqrt(2))).^2;
                        case 'pow3'
                            ctrl=(mean((pwr(h,4)))-3).^2;
                    end
                    if sum(ctrl)>table1(i)+table1(j)
                        %bad extrem indicated
                        rotated([i j])=1; %do not test the rotated signals anymore
                        W(:,[i j])=W(:,[i j])*rot2d;
                        repeat=1; %continue in iterating - the test of saddle points is positive
                        MaxIt=30;
                    end
                end
            end
        end
    end %if SaddleTest
end %while repeat
 
W1=W'*C^(-1/2);

%one-unit stage (for all sources)
for k=1:dim
    w=W(:,k);
    w=w/norm(w);
    wold=zeros(dim,1);
    NumIt=0;
    while (abs(w'*wold)<1-epsilon && NumIt<100 && ...
            (abs((W(:,k)/norm(W(:,k)))'*(w/norm(w)))>0.90))
        wold=w;
        u=w'*x;
        switch g
            case 'tanh'
                hyptan=tanh(u);
                w=(x*hyptan' - sum(1-hyptan.^2)'*w)/N;
            case 'gaus'
                gauss=u.*exp(-u.^2/2);
                dGauss=(1-u.^2).*exp(-u.^2/2);
                w=(x*gauss'-sum(dGauss)*w)/N;
            case 'biga'
                u=u';   
                m=mean(abs(u)); %estimate centers of distribution's 
                e=sqrt(1-m^2); %then their variance is..., because the overall variance is 1
                if e<=0.05
                    e=0.05; 
                    m=sqrt(1-e^2); 
                end %due to stability
                uplus=u+m; uminus=u-m;
                expplus=exp(-uplus.^2/2/e^2);
                expminus=exp(-uminus.^2/2/e^2);
                expb=exp(-(u.^2+m^2)/e^2);
                gg=-(uminus.*expminus + uplus.*expplus)./(expplus+expminus)/e^2;
                gprime=-(e^2*(expplus.^2+expminus.^2)+(2*e^2-4*m^2)*expb)./(expplus+expminus).^2/e^4;
                w=x*gg/N-mean(gprime)*w;                   
            case 'pow3'
                w=(x*(pwr(u,3))')/N-3*w;
        end
        w=w/norm(w);
        NumIt=NumIt+1;
    end
    if (abs((W(:,k)/norm(W(:,k)))'*(w/norm(w)))>0.95)
        W(:,k)=w;
    end
end

W2=W'*C^(-1/2);

function x=pwr(a,n)
x=a;
for i=2:n
    x=x.*a;
end