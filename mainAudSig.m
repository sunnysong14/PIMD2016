function mainAudSig(PMD, StrICA_lst, P_noise, Seed, OptVar_dir_loadData)
% Usage:
%   (a)Main function to systematically run Image Seperation Experiments. 
%   (b)This SINGLE file can also be used to PLOT, by either having the saved
%   results in the default directory in this function or assigning the Optimal 
%   input $OptVar_dir_loadData$. For example, one can put (ONLY) this
%   function to the result directory, set >> "OptVar_dir_loadData='';",
%   configurate INPUTS according to existing ICA results, and run this
%   function. Then one can obtain the .fig graphs. 
%   Please note, when there are no existing results in the pointing directory, 
%   this function will do ICA calculation if all required function
%   included; whereas if ONLY this function is using, errors will be
%   evoked saying "No defined function/variable".
%   It is worthwhile to note that to speed up PLOT, we can manipulate 
%   "if 1/0 ... end" in the save&outprint subsection.
% 
%   This file can run multiple (a) pMD values (b) noise levels $p$, 
%   (c) real-data cases, (d) sample sizes, for (e) ICA methods, 
%   (f) across multiple seeds.
% 
%   The results of each noise level, each real-data case, and each sample size 
%   are saved into .mat. The output Excel is divided for each real-data on
%   each PMD value but across all noise levels for all ICAs.
% 
% INPUTS: 
%   PMD -- A vector of pMD, i.e. the $p$ parameter of Mean Displacement.
%   StrICA_lst -- a STR list of names $strICA$.
%   P_noise -- Vector with positive constants defining the noise powers with value between 0
% and 1. In math, $pNoise=\sigma^2$.
%   OptVar_dir_loadData -- Optimal variable that can be omitted to replace $to_dir$ within
%   the code. It represents the directory where ".mat" results are saved.
%   By mannually inputing this variable, we can use this function to ONLY
%   PLOT, e.g. pNoise~SINR_LOSS.
% 
% OUTPUTS: Valid for the last $pMD$ with the last $pNoise$ ONLY.
%   Amari -- Matrix of amari errors. The amari errors of each ICA algorithm
%   is in row, and the errors for each random seed is in column.
%   SINR -- Matrix of $mean_sinr$. Each algorithm is in row and the the
%   errors for each random seed is in columns. Here, $mean_sinr$ is the
%   average sinr_ave across all observation signals.
%   SINRBound -- Matrix of upper bound of mean_SINR (across all signals), 
%   i.e. SINR-opt. The same with SINR: ICA algorithms are in rows, and the 
%   columns correspond to the seeds. 
%   Ave_amari -- Column vector of the mean amari for each ICA across all seeds.
%   Med_amari -- Column vector, median.
%   Std_amar -- Column vector, STD.
%   Ave_sinr, Med_sinr, Std_sinr -- the same as for Amari.
%   Ave_runtime, Std_runtime, Med_runtime -- run time counted in "tic-toc".
%   
%   %---------------------------------------------------------------%
%   The following parameters are only useful when length(Seed)=1. 
%   If there are multiple seeds, they are for the last seed.
%   %---------------------------------------------------------------%
%   S -- Spimdce
%   A -- initial mixing matrix
%   X_pn -- $p$-level noise corrupted X
%   X -- noise-free X
%   S_est -- estimated S
%   W_est -- estimated inverse of A, i.e. demixing matrix.
% 
% [NOTE] 
%   The INVALID Amari and SINR are handled here denoted by '-1'.

% [sample] Play audio:
% k=3;sig=S_est(k,:);Fs=16000;player=audioplayer(sig,Fs);play(player)

%% Prepare
% dir of data
if ~exist('OptVar_dir_loadData', 'var') %Default dir that ICA results are saved.
    to_dir = ['..', filesep, 'rslt', filesep];
else
    to_dir = OptVar_dir_loadData;
end 
if ~exist(to_dir, 'file'), mkdir(to_dir); end

% Debug: convert a single ICA to cell format
if ischar(StrICA_lst), StrICA_lst = cellstr(StrICA_lst); end

nSd = length(Seed);
nICA = length(StrICA_lst);
% screen
tag_verbose = 1;
shwIter = 10;
%% Main loop
for pmd = 1:length(PMD)
    pMD = PMD(pmd);
    if tag_verbose, fprintf('pMD: %0.2f:\n', pMD); end %screen output    
    
    % Each pNoise
    for pn = 1:length(P_noise)
        pNoise = P_noise(pn);
        if tag_verbose, fprintf('\tp_noise: %0.2f:\n', pNoise); end
        
        % name
        name = ['audSig_pN',strrep(num2str(pNoise),'.','_'),'_nSd',num2str(length(Seed)),'_ICA_'];
        for ii = 1 : length(StrICA_lst)
            strICA = StrICA_lst{ii};
            name = [name, strICA(1:3), '_'];
        end
        name = [name, 'pMD', strrep(num2str(pMD), '.', '_')];
        flnm = [to_dir, name];

        % init
        SINR = -ones(nICA, nSd);
        SINRBound = -ones(nICA, nSd);
        Error_innerC_realA = -ones(nICA, nSd);
        Error_innerC_estA = -ones(nICA, nSd);

        % Load or Calculate
        if exist([flnm, '.mat'], 'file')
            load(flnm);
            if tag_verbose, fprintf('\t\t Loading all ICAs, No Calculation ...\n'); end
        else
            %% Doing ALL ICAs across ALL seeds            
            % Run all ICAs seed-by-seed
            for sd = 1 : nSd
                seed = Seed(sd);

                % Obtain data
                [X_pn,~,~,A,Noise_cov] = obtainICA_signal(pNoise,seed);

                % Each ICA
                for ic = 1:nICA
                    strICA = StrICA_lst{ic};
                    % screen
                    if tag_verbose&&seed==1, fprintf('\t\t Calculating ICA: %s...\n', strICA); end
                    if ~mod(seed,shwIter) && tag_verbose, fprintf('\t\t\t%s:... seed: %d ...\n', strICA, seed); end

                    % [core]
                    [~,W_est,A_est,C] = mainICA_fun(strICA,pMD,X_pn,Noise_cov,A);

                    % check orthogonality of C with real A
                    switch lower(strICA)
                        case {'pegi', 'pimd'}
                            Error_innerC_realA(ic,sd) = compErrorInnerC(A,C);
                        otherwise  % invalid C
                            Error_innerC_realA(ic,sd) = -1;
                    end%END OF SWITCH

                    % check orthogonality of C with A_est
                    switch lower(strICA)
                        case {'pegi', 'pimd'}
                            Error_innerC_estA(ic,sd) = compErrorInnerC(A_est, C);
                        otherwise  % invalid C
                            Error_innerC_estA(ic,sd) = -1;
                    end%END OF SWITCH

                    % Calculate theoretical optimal demixing matrix for SINR-type
                    % things. It is used to calculate sinrbound_ave.
                    if min(svd(Noise_cov)) > 0
                        temp = (Noise_cov)^(-1/2) * A;
                        W_SINR_OPT = (temp'/(temp*temp'+eye(size(A))))*Noise_cov^(-1/2);
                    else % We assume noise_cov=0, and it is true in pimd experiment.
                        W_SINR_OPT = inv(A);
                    end

                    % Convergence, esp. fastICA-variants: 
                    % We enlarge the $max-iteration$ and $re-start$ for fastICA algorithms
                    % to enhance their convergence. It may help here since the real-data 
                    % kurtosis is often large enough. As still the inconvergence
                    % happens, we use "-1" to represent the unconvergence case.
                    % [Later] Add var $W_est_rm0row$. GIICA may not converge,
                    % which obtains only one or two rows of (e.g.) A and W.
                    W_est_rm0row = W_est();
                    W_est_rm0row(all(W_est_rm0row==0, 2), :) = [];
                    if size(W_est_rm0row,1) < size(A, 1)
                       sinr_ave = -1;
                    else
                        % Mean SINR across all estimated signals, i.e. all dimensions.
                        % it should be always valid.
                        % The following steps are to simplify calculation of SINR
                        % and should only change the demixing order.
                        W_est2 = unpermute(A, W_est);
                        sinr_ave = SINR_mean(A, W_est2, Noise_cov);
                    end
                    % sinrbound_ave exists anyway
                    sinrbound_ave = SINR_mean(A, W_SINR_OPT, Noise_cov);
                    SINR(ic, sd) = sinr_ave;
                    SINRBound(ic, sd) = sinrbound_ave;
                end%END OF FOR-ICA  
            end%END OF FOR-SEED

            % PF calulate Ave 
            Ave_sinr = -ones(size(SINR, 1), 1);
            Ave_sinrloss = -ones(size(SINR, 1), 1);
            % row-by-row
            for ii = 1 : size(SINR,1)
                sinr_arr = SINR(ii,:);
                % SINR
                Ave_sinr(ii) = mean(sinr_arr);
                % SINR-loss
                sinrBound = SINRBound(ii,:);
                sinrloss = sinrBound-sinr_arr;
                Ave_sinrloss(ii) = mean(sinrloss);
            end%END OF FOR-II
            % SINRBound
            Ave_sinrbound = mean(SINRBound,2);
            % error_innerC_realA
            Ave_error_innerC_realA = mean(Error_innerC_realA,2);
            Ave_error_innerC_estA = mean(Error_innerC_estA,2);
        end%END OF IF-Exist-load ELSE-ICA        
        %% [EXTRA Assianment]        
        % PLOT - Give pMD: pNoise~SINR_LOSS
        PlotPNoise_SinrLoss(:,pn) = Ave_sinrloss;
        
        % PLOT - Given pNoise: pMD~SINR_LOSS
        PlotPMD_SinrLoss_celPNoise{pn}(pmd,:) = Ave_sinrloss'; % matrix in cell
        %% SAVE & PRINT        
        % save .mat
        if 1
        save(flnm, 'StrICA_lst', 'pNoise', 'Seed', ...
            'SINR', 'SINRBound', 'Error_innerC_realA', 'Error_innerC_estA', ...
            'Ave_sinr','Ave_sinrloss', 'Ave_sinrbound', ...
            'Ave_error_innerC_realA', 'Ave_error_innerC_estA');
        end%END 1/0

        % excel printf
        nameExcel = ['audioSig','_pMD', strrep(num2str(pMD), '.', '_'),'.ods'];
        flnm_excel = [to_dir, nameExcel];
        hf = fopen(flnm_excel, 'a+');
        % info.
        fprintf(hf, [repmat('%s\t',1,8),'\n'], ...
            'nRun', 'pNoise', 'pMD',...
            'ICA','MeanSINR', 'SINRLoss', ...
            'ErC_realA','ErC_estA');
        % results:
        for ii = 1 : length(StrICA_lst)
            strICA = StrICA_lst{ii};
            fprintf(hf,'%d\t %f\t %d\t',length(Seed),pNoise,pMD);
            fprintf(hf,['%s\t',repmat('%f\t',1,5),'\r\n'], ...
                strICA, Ave_sinr(ii), Ave_sinrloss(ii), ...
                Ave_error_innerC_realA(ii), Ave_error_innerC_estA(ii));
            fprintf(hf,'\n');
        end
        fprintf(hf, '\r\n');
        fclose(hf);
    end%END OF p-NOISE
    %% PLOT pNoise~SINR_LOSS, given pMD
    if 1
    % Settings for pp-wt
    LineWidth = 10;
    MarkerSize = 20;
    MarkerFacecolor = [0.5,0.5,0.5];
    FontSize = 30;
    LegendLocation = 'northeast';   

    figure, 
    ID_legPltICA = 1; % Legend: record the ICA algo that is plotted.
    for k = 1:length(StrICA_lst)
        iplotICA = StrICA_lst{k};
        iplotICAsinrloss = PlotPNoise_SinrLoss(k,:);

        % Skip ICAs that are not plotted
        SkipICA_cel = {'ainv','sinropt','fasticatanh'};
        if ~isempty(strcmpi(SkipICA_cel, iplotICA))
            
            % line & marker for ICAs
            switch lower(iplotICA)
                case lower('PIMD')
                    LineMarker = 'diamond k-';
                    ppwt_iplotICA = 'PIMD';
                case lower('PEGI')
                    LineMarker = '> c-';
                    ppwt_iplotICA = 'PEGI+MD';
                case lower('GIICA')
                    LineMarker = 's m-';
                    ppwt_iplotICA = 'GIICA+MD';
                case lower('1fica')
                    LineMarker = 'x r-.';
                    ppwt_iplotICA = '1FICA+MD';
                case lower('fastica')
                    LineMarker = '*-.';
                    ppwt_iplotICA = 'FastICA+MD';
                case lower('fasticagauss') % robustICA
                    LineMarker = '^-.';
                    ppwt_iplotICA = 'rFastICA+MD';
                case lower('jade')
                    LineMarker = 'o:';
                    ppwt_iplotICA = 'JADE+MD';
                case lower('infomax')
                    LineMarker = '>-.';
                    ppwt_iplotICA = 'Informax+MD';
                otherwise
                    error(['Plot properties of ICAName "', iplotICA, '" have not been defined.']);
            end%END OF SWITCH
            
            % plot
            semilogy(P_noise,iplotICAsinrloss,LineMarker,...
                'LineWidth',LineWidth,'MarkerSize',MarkerSize,'MarkerFaceColor', MarkerFacecolor);

            % For legend
            legPlotICA{ID_legPltICA} = ppwt_iplotICA; ID_legPltICA = ID_legPltICA+1;
            hold on            
        else
        end        
    end%END OF FOR-plotICA

    % info
    L = legend(legPlotICA);
    L.Location = LegendLocation;
    xlabel('Noise Level')
    ylabel('SINR LOSS')
    title(strrep(['audSig', ', pMD=', num2str(pMD)], {'_'}, '-'));

    % set for pp-wr
    set(gca,'FontSize', FontSize);
    ylim([0, 2])
    end%END OF IF-PLOT: pNoise~pf    
end%END OF FOR-PMD
%% PLOT pMD~SINR_LOSS, given pNoise
if 1
% Settings for pp-wt
LineWidth = 10;
MarkerSize = 20;
MarkerFacecolor = [0.5,0.5,0.5];
FontSize = 30;
LegendLocation = 'northeast';

% One graph for each pNoise
for pn = 1 : length(P_noise)
    pNoise = P_noise(pn);
    PlotPMD_SinrLoss = PlotPMD_SinrLoss_celPNoise{pn}'; % Note transpose
    
    figure, 
    ID_legPltICA = 1; % Legend: record the ICA algo that is plotted.
    for k = 1 : length(StrICA_lst)
        iplotICA = StrICA_lst{k};
        iplotICAsinrloss = PlotPMD_SinrLoss(k,:);

        % Skip ICAs that are not plotted
        SkipICA_cel = {'ainv','sinropt','fasticatanh'};
        if ~isempty(strcmpi(SkipICA_cel, iplotICA))
            
            switch lower(iplotICA)
                case lower('PIMD')
                    LineMarker = 'diamond k-';
                    ppwt_iplotICA = 'PIMD';
                case lower('PEGI')
                    LineMarker = '> c-';
                    ppwt_iplotICA = 'PEGI+MD';
                case lower('GIICA')
                    LineMarker = 's m-';
                    ppwt_iplotICA = 'GIICA+MD';
                case lower('1fica')
                    LineMarker = 'x r-.';
                    ppwt_iplotICA = '1FICA+MD';
                case lower('fastica')
                    LineMarker = '*-.';
                    ppwt_iplotICA = 'FastICA+MD';
                case lower('fasticagauss') % robustICA
                    LineMarker = '^-.';
                    ppwt_iplotICA = 'rFastICA+MD';
                case lower('jade')
                    LineMarker = 'o:';
                    ppwt_iplotICA = 'JADE+MD';
                case lower('infomax')
                    LineMarker = '>-.';
                    ppwt_iplotICA = 'Informax+MD';
                otherwise
                    error(['Plot properties of ICAName "', iplotICA, '" have not been defined.']);
            end
            
            % plot
            semilogy(PMD,iplotICAsinrloss,LineMarker,...
                'LineWidth',LineWidth,'MarkerSize',MarkerSize,'MarkerFaceColor', MarkerFacecolor);

            % For legend
            legPlotICA{ID_legPltICA} = ppwt_iplotICA; ID_legPltICA = ID_legPltICA+1;
            hold on
        else
        end        
    end%END OF EACH-ICA

    % info
    L = legend(legPlotICA);
    L.Location = LegendLocation;
    xlabel('$p$-MD', 'Interpreter','latex')
    ylabel('SINR LOSS')
    title(strrep(['audSig, pNoise=', num2str(pNoise)], {'_'}, '-'));

    % set for pp-wr
    set(gca,'FontSize', FontSize);
    ylim([0, 2])    
end%END OF PN    
end%END OF IF-PLOT: pMD~pf
end%END OF FUNCTION
