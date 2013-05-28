function gi = getis(w,x)
% Getis statistic and global Moran's I for filtered variable
% -------------------------------------------------------------------------
% USAGE:
% gi = getis(w,x)
% mit:    w = weights matrix
%         x = vector with variable
% -------------------------------------------------------------------------
% OUTPUT:
% a structure variable:
%            gi.g = G_i statistic
%       gi.star_g = G_i* statistic
%          gi.e_g = E[G_i]
%     gi.star_e_g = E[ G_i*]
%           gi.y1 = Y1-parameter for variance of G_i
%      gi.star_y1 = Y1-parameter for variance of G_i*
%           gi.y2 = Y2-parameter for variance of G_i
%      gi.star_y2 = Y2-parameter for variance of G_i*
%          gi.s_g = standard deviation of G_i
%     gi.star_s_g = standard deviation of G_i*
%          gi.z_g = z-standardized G_i
%     gi.star_z_g = z-standardized G_i*
%          gi.x_f = filtered variable
%         gi.z_mi = z-standardized Moran's I 
%                   (of filtered variable)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Getis statistic and spatial filtering:
% -------------------------------------------------------------------------

n = length(w);                                       % number of observations

gi.g = (w*x)./(sum(x)-x);                            % G_i statistic

gi.e_g = (sum(w')/(n-1))';                           % Erwartungswert G_i

gi.y1 = (sum(x)-x)/(n-1);                            % Y1-parameter for variance of G_i       
gi.y2 = (sum(x.^2)-x.^2)./(n-1) - gi.y1.^2;          % Y2-parameter for variance of G_i
gi.s_g = sqrt(((sum(w').*(n - 1 - sum(w')))./...     % standard deviation of G_i
         ((n-1)^2*(n-2)))'.*(gi.y2./(gi.y1.^2)));
gi.z_g = (gi.g - gi.e_g)./gi.s_g;                    % z-standardized G_i

gi.x_f = (x.*gi.e_g)./gi.g;                          % filtered variable

gi.star_g = (w*x)./(sum(x));                                % G_i* statistic

gi.star_e_g = (sum(w')/(n-1))';                             % E[G_i*]

gi.star_y1 = (sum(x))/(n-1);                                % Y1-parameter for variance of G_i*     
gi.star_y2 = (sum(x.^2))./(n-1) - gi.star_y1.^2;            % Y2-parameter for variance of G_i*
gi.star_s_g = sqrt(((sum(w').*(n - 1 - sum(w')))./...       % standard deviation of G_i*
         ((n-1)^2*(n-2)))'.*(gi.star_y2./(gi.star_y1.^2)));
gi.star_z_g = (gi.star_g - gi.star_e_g)./gi.star_s_g;       % z-standardized G_i*


% -------------------------------------------------------------------------
% Global Moran's I for the filtered variable:
% -------------------------------------------------------------------------

y = gi.x_f;                             % endogenous variable
X(1:n,1) =1;                            % exogenous variable (a constant)

k = size(X);                            % dimension of X
k = k(:,2);                             % number of columns of X
df = n - k;                             % degrees of freedom
M = eye(n) - X*inv(X'*X)*X';            % residual maker
My = M*y;                               % residuals
denom = My'*My;                         % denominator of MI
mi = (My'*w*My)/denom;                  % Moran's I


K = M*0.5*(w+w')*M;                     % moments with trace-formulation:
ex_mi = (trace(K))/df;                          % expected value
var_mi = (2*(df*trace(K^2)-trace(K)^2))/ ...    % variance
         (df^2*(df+2));

gi.z_mi = (mi - ex_mi)/sqrt(var_mi);    % z-standardized Moran's I