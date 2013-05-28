% Example E
% -------------------------------------------------------------------------
% finds optimal delta for objfct.m, optimal weights matrix and
% filtered variable
% -------------------------------------------------------------------------
% OUTPUT:
% a structure variable
%            gi.g = G_i statistic
%          gi.e_g = Erwartungswert
%           gi.y1 = Y1-parameter for variance
%           gi.y1 = Y2-parameter for variance
%          gi.s_g = standard deviation
%          gi.z_g = z-standardized G_i
%          gi.x_f = filtered Variable
%         gi.z_mi = z-standardized Moran's I 
%                   (of filtered variable)
%
%               w = row-sum standardized weights matrix
%               u = unstandardized weights matrix
% 
%           delta = optimal delta
% -------------------------------------------------------------------------

clear                            
clc                              

% ---------------------------------------------------------------------
% Minimize objective function:
% ---------------------------------------------------------------------

[delta] = fminbnd(@objfct,1,300);

% ---------------------------------------------------------------------
% Optimal weights matrix:
% ---------------------------------------------------------------------

load('y.txt', '-ascii');                      % load variable
x = y(:,1);                                                  
d = load('distanz.txt', '-ascii');            % load distance matrix
[w,u] = distance2weight(d,delta);             % convert to
                                              % weights matrix

% ---------------------------------------------------------------------
% Optimal filtered variable:
% ---------------------------------------------------------------------

gi = getis(w,x);                              % filter variable
