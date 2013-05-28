function mp = moranplotmatrix(X,w)
% Moran scatterplot matrix
% -------------------------------------------------------------------------
% USAGE:
% mp = moranplotmatrix(X,w)
% mit:    X = matrix with variables
%         w = spatial weights matrix
% -------------------------------------------------------------------------
% OUTPUT:
% figure with Moran scatterplot matrix
% a structure variable
%            mp.y = matrix with variables
%                  (as deviation from mean)
%           mp.wy = matrix with spatially lagged variables
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Error checking:
% -------------------------------------------------------------------------

if nargin ~= 2
    error('Wrong number of arguments for moranplot');
end;

dimB = size(w);
if dimB(1) ~= dimB(2)
  error('Weights matrix is not quadratic'); 
end

if length(X) ~= dimB(1)
    error('Number of elements in variable do not correspond to number of elements in weights matrix'); 
end

% -------------------------------------------------------------------------
% Calculations:
% -------------------------------------------------------------------------

dim = size(X);                                  % dimension of Y
nVar = dim(2);                                  % number of variables

for i = 1:nVar
    mp.y(:,i) = X(:,i) - mean(X(:,i));          % variables (as deviation from mean)
    mp.wy(:,i) = w*mp.y(:,i);                   % spatially lagged variables
end

s = 1;                                          % coordinate matrix
for i = 1:nVar                                  % for individual 
    for j = 1:nVar                              % scatterplots
        coord(s,1) = j;
        coord(s,2) = i;
        s = s + 1;
    end
end

for i = 1:nVar*nVar
    r = coord(i,1);                             
    p = coord(i,2);                             
    subplot(nVar,nVar,i), ...
        plot(mp.y(:,r),mp.wy(:,p), ...
        '.', ...                                
        'MarkerEdgeColor','r', ...             
        'MarkerSize',5);                        
    lsline;                                     
    xlabel(['y_',num2str(r)]);                  
    ylabel(['wy_',num2str(p)]) ;               
    grid on;                                  
end

