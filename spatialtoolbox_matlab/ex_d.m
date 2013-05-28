% Example D
% -------------------------------------------------------------------------
% Demonstrates usage of  moranplotmatrix.m
% -------------------------------------------------------------------------

clear                                           
clc                                             


y = load('y.txt', '-ascii');                    % load variable
y = y(:,1:2);

w = load('w.txt', '-ascii');                    % load weights matrix

mp = moranplotmatrix(y,w);                     

