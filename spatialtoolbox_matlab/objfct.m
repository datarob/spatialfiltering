function f = objfct(delta)
% Objective function for optimal delta (see ex_e.m)
% -------------------------------------------------------------------------

load('y.txt', '-ascii');                       % load variable
x = y(:,1);                                                
d = load('distanz.txt', '-ascii');             % load distance matrix

[w,u] = distance2weight(d,delta);              % convert to
                                               % weights matrix
gi = getis(w,x);                               % filter variable

f = abs(gi.z_mi);                              % absolut z-value of
                                               % global Moran's I