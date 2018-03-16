% This function is from fastICA, which would be treated as a general
% function to remove the mean of matrix X. 
% 
%  o The setting of "vectors" is consistent with my habit, that is row
%  corresponds to each dimension, and columns contain samplings.
% 
% Liyan 05-05-2016
%%
function [newVectors, meanValue] = remmean(vectors)
%REMMEAN - remove the mean from vectors
%
% [newVectors, meanValue] = remmean(vectors);
%
% Removes the mean of row vectors.
% Returns the new vectors and the mean.
%
% This function is needed by FASTICA and FASTICAG

% @(#)$Id: remmean.m,v 1.2 2003/04/05 14:23:58 jarmo Exp $

newVectors = zeros (size (vectors));
meanValue = mean (vectors')';
newVectors = vectors - meanValue * ones (1,size (vectors, 2));
