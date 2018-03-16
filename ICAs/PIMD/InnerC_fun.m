function C = InnerC_fun(X)
% Usage:
%   Specific function for OUR noisy ICA model. It aims to form our matrix of 
%   inner product. Our procedure consists of the following steps:
%    - (1) M1 = Q_x(I);
%    - (2) M2 = inv(M1);
%    - (3) M3 = Q_x{M_2};
%    - (4) C = inv(M3).
% 
%   Based on the Theorem 4.3 of [VossLine_COLT2013], $M_2 = Q_x(M_1) = ADA^T$,
%   where $D_{ii} = 1/||A_i||^2$ is positive. It is arrived that 
%   $C = (A^T)^{-1} D' A^{-1}$, where $D'_{ii} = ||A_i||^2$ is positive.
% 
% Input:
%   X -- centered X. Row vectors are dimension, and column vectors are
%   samplings.
% 
% Output:
%   C -- OUR inner product matrix related to a Basis of R^N space, under
%   which {A1, A2, ..., An} is an orthogonal basis where Ai is the i-th
%   column of (unwhitened) mixing matrix A.
% 
% Note:
%   We use Moore-Penrose pseudo-inverse instead of inverse.   
%   If A is square and not singular, then pinv(A) is an expensive way to 
%   compute inv(A); If A is not square, or is square and singular, then 
%   inv(A) does not exist. In these cases, pinv(A) has some of, but not all,
%   the properties of inv(A).
% 
% Reference: 
%   @inproceedings{VossLine_COLT2013,
%     author 		= 	{Belkin, M. and Rademacher, L. and Voss, J.},
%     title 		= 	{Blind Signal Separation in the Presence of Gaussian Noise},
%     booktitle 	= 	{Conference on Learning Theory (COLT)},
%     pages 		= 	{270-287},
%     year          = 	{2013}}
% 
% Liyan for NIPS16 05-04-2016

%%
Cum4_tsor = Cum4_fun(X);  % require centered X inherently
m = size(X, 1);  % All dimensions of $Cum4_tsor$ should be equal.

% step 1. compute Qx(I)
II = eye(m);
M1 = QxM_fun(Cum4_tsor, II);

% step 2. inverse M1
M2 = pinv(M1);
% M2 = inv(M1);

% step 3. compute Qx(M2)
M3 = QxM_fun(Cum4_tsor, M2);

% step 4. set $C = inv(M3)$
C = pinv(M3);
% C = inv(M3);

end  % END OF FUNCTION
