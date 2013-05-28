function [w,u] = distance2weight(d,delta)
% Convert a distance matrix to a weights matrix
% -------------------------------------------------------------------------
% USAGE:
% [w,u] = distance2weight(d,delta)
% with:    d = distance matrix
%     delta = delta-parameter for distance decay function
% -------------------------------------------------------------------------
% OUTPUT:
%            w = row-sum standardized weights matrix
%            u = unstandardized weights matrix
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Converting distance matrix to weights matrix:
% -------------------------------------------------------------------------

beta = 1/delta/1000;            % convert to kilometers
w = exp(-d*beta);               % distance decay function

n = length(w);                  % number of observations
w = w-eye(n);                   % zeros in main diagonal
u = w;                           % unstandardisierte Gewichtungsmatrix

% -------------------------------------------------------------------------
% Row-sum standardization:
% -------------------------------------------------------------------------

m = ndims(w);                  % number of array dimensions
nterm = sum(w,m);              % row-sum standardization
nterm = repmat(nterm,[ones(1,m-1) size(w,m)]);
nterm = nterm + (nterm==0);    % protect against zeros in division
w = w ./ nterm;                % row-sum standardized matrix
