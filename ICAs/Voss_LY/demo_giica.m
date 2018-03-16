% - Add Amari index part for each preprocessing method.
%
% Liyan for NIPS16 on 04-11-2016
%  last updated: 04-11-2016
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%function demo_giica()
% function demo_giica()
% This is a demo script for the GIICA routine.  It produces data (both the
% underlying clean set of data (uniform vs. uniform) and a perturbed
% noisy data that includes the additive Gaussian noise).  Then, the noisy
% data is run through GI-ICA after various preprocessing ideas including:
%   - k4 pseudo inner product -- This is not a preprocessing step, but
%       instead uses fourth cumulant information to set up an "inner
%       product space" (which may have negative values of < x, x >) in
%       which the cumulant based gradient iteration ICA can be run.  This
%       tends to produce more accurate results than
%       quasi-orthogonalization, and avoids some numerical issues.
%   - quasi-orthogonalization -- Fourth cumulant based preprocessing step
%       which approximately orthogonalizes the underlying independent
%       components.
%   - whitening -- Preprocesses the data using covariance information in
%       order to make it have identity covariance.  In the noise free
%       setting (no additive Gaussian noise), this is valid and preferred
%       for orthogonalizing the latent source signals.  However, it is not
%       invaraint to Gaussian noise.
% The demixing are then displayed with the demixing matrix applied both to
% the noisy data used by the ICA algorithm and to the underlying noise-free
% data.
%
% In practice, the user would only have access to the noisy data.

    % Toggle to 1 in order to repeat the previous run.  Note that this
    % saves to the same file as used by the other demo script.
    UsePriorRun = 0;

    if UsePriorRun
        load('vars')
        rng(status);
    else
        status = rng;
        rng(status);
    end
    %%save('vars')

    samples = 50000; %5000;
    NoiseMagnitude = 0.5;

    S = (rand(2, samples))*2*sqrt(3);
    A = [1 1
        0 1];
    Sigma = [1 1
        1 2] / 2.6180;
    X = A*S;
    X2 = X + mvnrnd(zeros(2, samples)', NoiseMagnitude*Sigma)';

    figure(1);
    plot(X2(1, :), X2(2, :), 'g+', X(1, :), X(2, :), 'k.');
    title('Raw Data');
    legend('noisy data', 'noiseless data');
    axis('equal');

    %% Test the GI-ICA interface function
    [S2, A_inv, b] = GIICA(X2, 'whiten', 'verbose', 3, 'SINR variant', 0);
     S = A_inv*(X - repmat(b, 1, samples));
    
    % No-1 normal whitening
    % Amari error: LY
    P = A_inv * A ;                       % permutation matrix
    id_ = 1; 
	Amari(id_) = Amari_index_ISA(P, [1,1], 'uniform', 2);  % traditional Amari-index 
     
    figure(2);
    plot(S2(1, :), S2(2, :), 'g+', S(1, :), S(2, :), 'k.');
    title({'demixed data (standard whitening & GI-ICA)',...
        'Demixing matrix determined using only the noisy data.'});
    legend('noisy data', 'noiseless data');
    axis('equal');

    %% Test the GI-ICA interface function
    [S2, A_inv, b] = GIICA(X2, 'quasi-orthogonalize',  'SINR variant', 0, 'verbose', 3);
     S = A_inv*(X - repmat(b, 1, samples));

    % No-2 GI-ICA, NIPS13
    % Amari error: LY
    P = A_inv * A ;                       % permutation matrix
    id_ = id_ + 1; 
	Amari(id_) = Amari_index_ISA(P, [1,1], 'uniform', 2);  % traditional Amari-index  
     
    figure(3);
    plot(S2(1, :), S2(2, :), 'g+', S(1, :), S(2, :), 'k.');
    title({'demixed data (quasi-orthogonalization & GI-ICA)',...
        'Demixing matrix determined using only the noisy data.'});
    legend('noisy data', 'noiseless data');
    axis('equal');

    %% Test the PEGI (pseudo-Euclidean gradient iteration) interface (4th cumulant)
    [S2, A_inv, b] = GIICA(X2, 'pseudo-Euclidean IP',  'SINR variant', 0, 'verbose', 3);
     S = A_inv*(X - repmat(b, 1, samples));

    % No-3 PEGI+normal
    % Amari error: LY
    P = A_inv * A ;                       % permutation matrix
    id_ = id_ + 1; 
	Amari(id_) = Amari_index_ISA(P, [1,1], 'uniform', 2);  % traditional Amari-index  
    
    figure(4);
    plot(S2(1, :), S2(2, :), 'g+', S(1, :), S(2, :), 'k.');
    title({'demixed data (k4 pseudo-Euclidean gradient iteration)',...
        'Demixing matrix determined using only the noisy data.'});
    legend('noisy data', 'noiseless data');
    axis('equal');

    %% Test the PEGI interface with SINR demixing (4th cumulant)
    [S2, A_inv, b] = GIICA(X2, 'pseudo-Euclidean IP', 'SINR variant', 1, 'verbose', 3);
     S = A_inv*(X - repmat(b, 1, samples));

    % No-4 PEGI+SINR
    % Amari error: LY
    P = A_inv * A ;                       % permutation matrix
    id_ = id_ + 1; 
	Amari(id_) = Amari_index_ISA(P, [1,1], 'uniform', 2);  % traditional Amari-index  
     
    figure(5);
    plot(S2(1, :), S2(2, :), 'g+', S(1, :), S(2, :), 'k.');
    title({'SINR optimal demixed data (k4 pseudo-Euclidean gradient iteration)',...
        'Demixing matrix determined using only the noisy data.'});
    legend('noisy data', 'noiseless data');
    axis('equal');

    %% Test the PEGI interface (3rd cumulant)
    X = A * exprnd(ones(2, samples));
    X2 = X + mvnrnd(zeros(2, samples)', NoiseMagnitude*Sigma)';

    [S2, A_inv, b] = GIICA(X2, 'pseudo-Euclidean IP', 'contrast', 'k3', 'SINR variant', 0, 'verbose', 3);
    S = A_inv*(X - repmat(b, 1, samples));

    % No-5
    % Amari error: LY
    P = A_inv * A ;                       % permutation matrix
    id_ = id_ + 1; 
	Amari(id_) = Amari_index_ISA(P, [1,1], 'uniform', 2);  % traditional Amari-index  
    Amari = Amari(:);
    
    figure(6);
    
    plot(S2(1, :), S2(2, :), 'g+', S(1, :), S(2, :), 'k.');
    title({'demixed data (k3 pseudo-Euclidean gradient iteration)',...
        'Demixing matrix determined using only the noisy data.', ...
        'Underlying sources follow the exponential distribution'});
    legend('noisy data', 'noiseless data');
    axis('equal');

% % end