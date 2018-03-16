% [08/11/2016] for PIMD'acml16
%   o updates original function with $pMD$. To change the codes as less as possible, 
% I use Global variable.
%  o add choice of preprocessing with ours as follow:
% "strcmp(preprocessing, 'real-euclidean ip'))
% 
% Outputs:
%   S: estimated sources, i.e. estimated ICs
%   W: estimated inverse of mixing matrix A
%   b: The mean of X which is subtracted during centering.
%   A: estimated mixing matrix, which follows W*A = I.
% Usage summary:
%   giica [NIPS13]: [S_est, W_est] = GIICA(X, 'quasi-orthogonalize',  'SINR variant', 0, 'verbose', 0);

%%
function [ S, W, b, A, C] = GIICA_LY(X, preprocessing, varargin)
%function [ S, W, b, A ] = GIICA( X, preprocessing, varargin )
%   All users should read at least to the dash line.
%   Required Inputs:
%   X:  The mixed data.  Each column of X represents a sample of data to be
%       demixed via ICA.  The ICA model is X = A*S (where S contains the
%       latent samples with independent coordinates, and A is the mixing
%       matrix).
%   preprocessing: -- designates the choice of method used in the first
%       step of the 2-step ICA algorithm.
%       Values:
%       'pseudo-Euclidean IP' -- Does not perform an explicit
%           preprocessing, but instead uses third or fourth order cumulant
%           methods (matching the choice of contrast) to generate a
%           psuedo inner product space in which to perform gradient
%           iteration ICA algorithm.  For details, see the NIPS 2015 paper
%           "A Pseudo-Euclidean Iteration for Optimal Recovery in Noisy
%           ICA" by Voss, Belkin, and Rademacher.
%       'quasi-orthogonalize' -- Uses fourth order cumulant 
%           methods to orthogonalize the latent signals.  This method 
%           is robust to arbitrary, additive Gaussian noise given
%           sufficient samples.  The 'pseudo-Euclidean IP' preprocessing is
%           the recommended preprocessing for the noisy setting.  This
%           algorithm is based on the 2013 NIPS paper "Fast Algorithms for
%           Gaussian Noise Invariant Independent Component Analysis" by
%           Voss, Rademacher, and Belkin for details.
%       'whiten' -- Standard whitening more traditionally used in ICA.
%           This technique should be used on noise-free data, as it
%           requires fewer samples than the noise-robust preprocessing
%           techniques to achieve good sample estimates.
%       'none' -- This forgos the preprocessing step.  Note that the 
%           ICA algorithms implemented are only valid for data which
%           is centered and quasi-orthogonalized (which is a weaker
%           condition than whitening on noise free data).  Only use this
%           option if the data already meets these constraints.
%
%   Return Values:
%   S:  The demixed data (S = W * (X - b)).
%   W:  The demixing matrix used.  By default, this attempts to construct
%       the SINR optimal demixing of the source signals (demixing using
%       A'*cov(X)).  This coincides with inv(A) only in the whitening case.
%       Set the SINR Variant flag to 0 to get demixing with inv(A) instead.
%   b:  The mean of X which is subtracted during centering.
%   A:  The approximate mixing matrix recovered under the ICA assumptions.
%
%   Optional Inputs:
%   GIICA should be called in the form:
%       GIICA(X, preprocessing, <op string 1>, <op 1 value>, <op string 2>, <op 2 value>, ... )
%   Where the elipses encompass that any number of options can be used.
%   All options are given default values.  Some example calls are:
%       GIICA(X, 'quasi-orthogonalize');
%       GIICA(X, 'whiten');
%       GIICA(X, 'quasi-orthogonalize', 'contrast', 'k3');
%       GIICA(X, 'pseudo-Euclidean IP', 'SINR variant', 1, contrast', 'k4');
%   The first format is for noisy data (with an unknown, 
%   additive Gaussian noise), and the second format is recommended when
%   the data is believed to be noise-free, or has fewer samples than
%   required to run the 'quasi-orthogonalize' version accurately.  The
%   third format uses an optional input to switch to a third-cumulant
%   skew-based implementation of ICA.  The third cumulant is useful when
%   the latent source signals are not symmetric.
%--------------------------------------------------------------------------
%   The possible optional inputs (and default values) are specified below:
%   Options:
%   'contrast' -- designates the function used during the second step of
%       the ICA algorithm.  This function will be used in the
%       gradient-iteration fashion to recover columns of the mixing matrix
%       under orthogonality assumptions.
%       Values:
%       'k3' -- Third cumulant, or the skew.
%       'k4' (default) -- Fourth cumulant.
%       'rk3' -- Welling's robust third cumulant.  Not recommended when
%           there is an additive Gaussian noise.  Robust in the sense of
%           outliers, but requires more samples than the traditional third
%           cumulant.
%       'rk4' -- Welling's robust fourth cumulant.  Not recommended when
%           there is an additive Gaussian noise.  Robust in the sense of
%           outliers, but requires more samples than the traditional fourth
%           cumulant.
%   'robustness factor' -- Used to specify the robustness constant alpha 
%       used by Welling's robust cumulants.
%       Value:  Provide a decimal value.  Default is 2.  Valid values range
%           from 1 (mimicks traditional cumulants) to infinity.
%       example usage:  GIICA(X, 'whiten', 'contrast', 'rk4', 'robustness factor', 1.5)
%   'SINR variant' -- Used to specify wether we are performing standard
%       ICA or the offline maximizer of signal to interference plus
%       noise ration version of ICA.  This distinction is only relevant
%       in the presence of an additive (assumed Gaussian) noise.
%       This option has no effect when using whitening as the
%       preprocessing.  The valid preprocessing choices for this to give
%       good results are 'pseudo-Euclidean IP' and 'quasi-orthogonalize'.
%       When performing whitening, this has no effect.
%       Value:  Boolean value of 0 or 1.  1 (default) indicates that we are
%           performing the SINR optimal version of ICA.  In this case,
%           the demixing matrix W is not equal to inv(A) in general.
%       example usage:  GIICA(X, 'pseudo-Euclidean IP', 'SINR variant', 1)
%   'max iterations' -- Specifies the cap on the number of iterations which
%       can be used to find a single column of A.
%       Value:  Integer value.  Default is 1000.
%       Example usage:  GIICA(X, 'quasi-orthogonalize', 'max iterations', 100)
%   'orthogonal deflation' -- Determines in the second step of ICA (the
%       deflationary recovery of the demixing matrix under an orthogonality
%       constraint) whether orthogonality is enforced between the column of
%       A currently being recovered and all previously recovered columns.
%       The algorithms employed are in theory self-orthogonalizing, but
%       this is not necessarily true on sample data or when the ICA
%       assumptions do not fully hold.
%       Value: 
%       'false' -- the orthogonality constraint is relaxed after
%           convergence is achieved in the orthogonal subspace.  The
%           maximum iterations bounds the sum of orthogonality enforced and
%           orthogonality relaxed iterations.  If convergence is never
%           achieved in the orthogonal subspace, then orthogonality remains
%           enforced despite this flag.
%       'true' (default)
%   'precision'
%       Value:  A decimal value.  During the deflationary gradient
%           gradient iteration step, this value will be used by the
%           stopping criterion.  In particular, if 2 subsequent estimates
%           for a column of the mixing matrix differ (in terms of cosine)
%           by less than the specified precision, the routine finishes.
%           This value must be less than 1.  The default is 0.0001
%       example call:  GIICA(X, 'quasi-orthogonalize', 'precision', 1e-6)
%    'verbose' -- An integer value specifying the level of verbosity.
%       Value:
%       0  -- Displays runtime errors and MatLab generated Warnings only
%           (not recommended).
%       1  -- Displays runtime errors and warnings.
%       2 (default) -- Displays runtime errors, warnings, and high level 
%           progress indicators.
%       3  -- Displays runtime errors, warnings, and gives lower level
%           progress indicators.
%       Example usage:  GIICA(X, 'whiten', 'verbose', 1)
%
%   NOTE:  Where numerical values are expected for the input, we make no
%   guarantees that the error checking that the input is valid is fully
%   carried out.
%
%   NOTE:  This code is largely based on the works:
%   Fast Algorithms for Gaussian Noise Invariant Independent Component Analysis
%   by:  James Voss, Luis Rademacher, and Mikhail Belkin
%   (NIPS 2013)
%
%   A Pseudo-Euclidean Iteration for Optimal Recovery in Noisy ICA
%   by:  James Voss, Mikhail Belkin, and Luis Rademacher
%   http://arxiv.org/abs/1502.04148
    
    if (mod(nargin, 2) ~= 0) || (nargin < 2)
        s = sprintf('%s\n%s\n%s', ...
            'Error:  Invalid number of arguments', ...
            'BASIC USAGE:  [S, A_inv, b] = GIICA(X, preprocessing)', ...
            'ADVANCED USAGE:  [S, A_inv, b] = GIICA(X, preprocessing, <op string 1>, <op 1 val>, <op string 2>, <op 2 val>' );
        error(s);
    end
        
    %% Set default parameter values
    alpha = 2; % robustness factor
    contrast = 'k4';
    epsilon = 1e-4; % precision
    orthogonality = int32(1); % orthogonal deflation flag
    maxIterations = int32(1000);
    verbosity = int32(2);
    SINROptFlag = int32(1);

    %% Error checking for preprocessing choice
    preprocessing = lower(preprocessing);
    if ~( strcmp(preprocessing, 'none') || ...
            strcmp(preprocessing, 'quasi-orthogonalize') || ...
            strcmp(preprocessing, 'whiten') || ...
            strcmp(preprocessing, 'pseudo-euclidean ip') || ...
            strcmp(preprocessing, 'real-euclidean ip'))  % LY: add the last [05-06-16]
        error('Usage:  Invalid preprocessing choice:  %s', preprocessing);
    end
    
    %% Switch parameter values based on user input arguments
    numUserOps = (nargin - 1) / 2;
    for i = 0:numUserOps-1
        option = varargin(2*i + 1);
        option = option{1};
        opVal = varargin(2*i + 2);
        opVal = opVal{1};
        switch option
            case 'contrast'
                opVal = lower(opVal);
                if strcmp(opVal, 'k3') || strcmp(opVal, 'k4') || ...
                        strcmp(opVal, 'rk3') || strcmp(opVal, 'rk4')
                    contrast = opVal;
                else
                    error('Usage:  Invalid contrast function choice:  %s', opVal);
                end
            case 'robustness factor'
                if ~isnumeric(opVal)
                    error('Usage:  robustness factor must have a numerical type.');
                end
                if opVal >= 1
                    alpha = opVal;
                else
                    error('Usage:  robustness factor must be at least 1');
                end
            case 'max iterations'
                if ~isnumeric(opVal)
                    error('Usage:  Max Iterations must have a numerical type.');
                end
                opVal = int32(floor(opVal));
                if opVal <= 0
                    error('Usage:  Max Iterations must be strictly positive.');
                else
                    maxIterations = opVal;
                end
            case 'orthogonal deflation'
                opVal = lower(opVal);
                if strcmp(opVal, 'true')
                    orthogonality = 1;
                elseif strcmp(opVal, 'false')
                    orthogonality = 0;
                else
                    error('Usage:  Invalid choice for orthogonality constraint:  %s', opVal );
                end                
            case 'precision'
                if ~isnumeric(opVal)
                    error('Usage:  precision value must have a numerical type.');
                end
                if opVal < 1
                    epsilon = opVal;
                else
                    error('Usage:  Invalid precision value (must be less than 1):  %f', opVal);
                end
            case 'SINR variant'
                if opVal ~= 0 && opVal ~= 1
                    error(['Usage:  SINR variant must be an integer either ' ...
                           '0 or 1']);
                end
                SINROptFlag = opVal;                    
            case 'verbose'
                if ~isnumeric(opVal)
                    error('Usage:  verbosity value must have a numerical type.');
                end
                opVal = int32(floor(opVal));
                if opVal < 0 || opVal > 3
                    error('Usage:  Invalid parameter for ''verbose'' options');
                end
                verbosity = opVal;
            otherwise
                error('Usage:  Invalid option choice:  %s', option);
        end
    end

    %% Run the ICA algorithm given user specificaitons.
    [S, W, A, b, totSteps, C] = ...
        ICA_Implementation_LY( X, contrast, preprocessing, orthogonality, ...
                            epsilon, maxIterations, alpha, ...
                            verbosity, SINROptFlag);
end

