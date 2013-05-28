% Example F
% -------------------------------------------------------------------------
% Calculate weights matrix and spatially filterd variable
% -------------------------------------------------------------------------
% OUTPUT:
% a structure variable
%            gi.g = G_i statistic
%          gi.e_g = expected value
%           gi.y1 = Y1-parameter for variance
%           gi.y1 = Y2-paramteter for variance
%          gi.s_g = standard deviation
%          gi.z_g = z-standardized G_i
%          gi.x_f = filtered variable
%         gi.z_mi = z-standardized Moran's I 
%                   (of filtered variable)
%
%               w = rowsum-standardized weights matrix
%               u = unstandardized weights matrix
% 
% -------------------------------------------------------------------------

clear                                         
clc                                           

load('y.txt', '-ascii');                      % load variable
x = y(:,1);                                                
d = load('distanz.txt', '-ascii');            % load distance matrix

[w,u] = distance2weight(d,138.1611);          % convert to 
                                              % weights matrix

gi = getis(w,x);                              % filter variable