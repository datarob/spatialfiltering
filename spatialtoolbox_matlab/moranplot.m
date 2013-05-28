function[mp,mp_regstats] = moranplot(y,w)
% Moran scatterplot (with outlier detection)
% -------------------------------------------------------------------------
% USAGE:
% [mp,mp_regstats] = moranplot(y,w);
% with:    y = vector with variable
%         w = spatial weights matrix
% -------------------------------------------------------------------------
% OUTPUT:
% Figure with Moran scatterplot (linear and LOWESS),
% studentized residuals, Cook's distances
% a structure variable
%            mp.y = vector with variable
%                  (as deviation from mean)
%           mp.wy = vector with spatially lagged variable
%     mp_regstats = regression diagnostics
%                   (see also 'Statistics Toolbox' - regstats)   
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

if length(y) ~= dimB(1)
    error('Number of elements in variable do not correspond to number of elements in weights matrix'); 
end

% -------------------------------------------------------------------------
% Calculations:
% -------------------------------------------------------------------------

mp.y = y - mean(y);                          % variable ((as deviation 
                                             % from mean)

mp.wy = w*mp.y;                              % spatially lagged variable           

mp_regstats = regstats(mp.wy,mp.y, ...       % regression diagnostics
              'linear','all');      

sp1 = subplot(2,2,1);
plot (mp.y,mp.wy, ...                        % Moran scatterplot (linear)
      '.', ...                               % points
      'MarkerEdgeColor','r', ...             % color of points
      'MarkerSize',5)                        % size of points
lsline                                       % linear regression line
title('Moran scatterplot (linear)')         
xlabel(['y, ',...                           
        ' MI = ', ...
        num2str(mp_regstats.beta(2))])       % slope of regression line
ylabel('wy')                                 
grid on                                      


yy = smooth(mp.y,mp.wy,0.4,'lowess');        % LOWESS smoother, span = 40%
[xx,ind] = sort(mp.y);                       % sort data for figure

sp2 = subplot(2,2,2);
plot(mp.y,mp.wy,'r.',xx,yy(ind),'b-')        % Moran scatterplot (LOWESS)  
title('Moran scatterplot (LOWESS)')          
xlabel('y')                                  
ylabel('wy')                                 
grid on                                      
                
sp3 = subplot(2,2,3);                        
plot(mp_regstats.studres)                    % studentized residuals
title('Studentized residuals')             
xlabel('Observation')                        
set(sp3,'YGrid','on')                        

sp4 = subplot(2,2,4);
plot(mp_regstats.cookd)                      % Cook's distance
title('Cook''s distance')                   
xlabel('Observation')                        
set(sp4,'YGrid','on')                       

