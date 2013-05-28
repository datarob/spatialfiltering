% Example B
% -------------------------------------------------------------------------
% Demonstrates usage of moranplot.m
% -------------------------------------------------------------------------

clear                                           
clc                                             

y = load('y.txt', '-ascii');                    % load variable
y = y(:,1);

w = load('w.txt', '-ascii');                    % load weights matrix

[mp,mp_regstats] = moranplot(y,w);             
