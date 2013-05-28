% Example A
% -------------------------------------------------------------------------
% Demonstrates usage of moransad.m
% -------------------------------------------------------------------------

clear                            
clc                             

y = load('y.txt', '-ascii');     % load endogenous variable                                
y = y(:,1);
nObs = length(y);                % number of observations

X(1:nObs,1) =1;                  % exogenous variable (a constant)
            
B = load('u.txt', '-ascii');     % load weights matrix

code = 'W';                      % coding scheme
GlobalMI = 'gl';                 % global (a. local MI)
Sad = 'j';                       % saddlepoint approximation

mi  = moransad(y,X,B,code,GlobalMI,Sad); 


