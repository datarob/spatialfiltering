function mi = moransad(y,X,B,code,GlobalMI,Sad);
% Global and local Moran's I (incl. saddle point approximation)
% -------------------------------------------------------------------------
% USAGE:
% mi = moransad(y,X,B,code,GlobalMI,Sad);
% with:    y = vector with endogenous variable
%         X = vector/matrix with exogenous variable(s)
%         B = spatial weights matrix
%      code = encoding method: 
%             c...globally standardized, 
%			  w...row-sum standardized,
%			  s...variance stabilizing 
%  GlobalMI = g... only global MI,
%             gl... global and local MI
%       Sad = j... perform saddle point approximation
%             n... don't perform saddle point approximation
% -------------------------------------------------------------------------
% OUTPUT:
% a structure variable
%            mi.My = residuals
%            mi.mi = Moran's I (first entry global, others local)
%   mi.moments_exp = expected value
%   mi.moments_var = variance
%   mi.moments_skw = skewness
%   mi.moments_kur = kurtosis
%          mi.z_mi = z-standardized MI
%       mi.prob_nv = probability with normal distribution
%           mi.sad = saddle point
%         mi.sad_r = r - parameter for distribution
%         mi.sad_u = u - parameter for distribution
%      mi.prob_sad = probability for saddlepoint approximation
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Error checking:
% -------------------------------------------------------------------------

if nargin ~= 6
    error('Wrong number of arguments for moransad');
end;

dimB = size(B);
if dimB(1) ~= dimB(2)
  error('Weights matrix is not quadratic');  
end

if length(y) ~= length(X)
    error('Number of elements in y do not correspond to X'); 
end

if length(y) ~= dimB(1)
    error('Number of elements in y does not correspond to spatial weights matrix'); 
end

warning off MATLAB:divideByZero         

tic;                                    % start timer

% -------------------------------------------------------------------------
% Create variables:
% -------------------------------------------------------------------------

nObs = length(y);                       % number of observations
B = sparse(B);                          % convert to sparse matrix
B = (B'+B)/2;                           % ensure symmetry

X = sparse(X);                          % convert to sparse matrix

if GlobalMI == 'g'                      % switch:
    nLoop = 1;                          % Global MI
else    
    nLoop = nObs + 1;                   % Global and local MI              
end    


k = size(X);                            % dimension of X
k = k(:,2);                             % number of columns of X
df = nObs - k;                          % degrees of freedom
M = eye(nObs) - X*inv(X'*X)*X';         % residual maker
mi.My = M*y;                            % residuals
denom = mi.My'*mi.My;                   % denominator in MI

% -------------------------------------------------------------------------
% Loop over global and local Moran's I:
% -------------------------------------------------------------------------

for i = 1:nLoop
    
    % ---------------------------------------------------------------------  
    % Encode weights matrix for global Moran's I:
    % ---------------------------------------------------------------------
    
    if i == 1
        if code == 'S'                      % s...variance stabilizing
            D = sparse(diag(sqrt(sum((B.*B)')')));
            D = inv(D);
            Vi = D*B;
            a = sum(sum(Vi')');
            Vi = nObs/a*Vi;
            Vi = 0.5*(Vi + Vi');            % ensure symmetry
        end
        if code == 'W'                      % w...row-sum standardized
            D = sparse(diag(sum(B')'));
            D = inv(D);
            Vi = D*B;
            Vi = 0.5*(Vi + Vi');            % ensure symmetry
        end
        if code == 'C'                      % c...globally standardized
            a = sum(sum(B')');
            Vi = nObs/a*B;
            Vi = 0.5*(Vi + Vi');            % ensure symmetry
        end  
        
    % ---------------------------------------------------------------------    
    % Encode weights matrix for local Moran's I:
    % ---------------------------------------------------------------------    
        
    else
        
        clc
        fprintf('Iteration %g of %g',i,nLoop)    
        
        Vi = sparse(1,nObs,0,nObs,nObs);    % star-shaped 
        Vi((i-1),:) = B((i-1),:);           % weights matrix
        Vi(:,(i-1)) = B(:,(i-1));
        if code == 'S'                      % s...variance stabilizing
            Vi = (nObs^2)*Vi*D((i-1),(i-1))/(2*a);
        end
        if code == 'W'                      % w...row-sum standardized
            Vi = nObs*Vi*D((i-1),(i-1))/2;
        end
        if code == 'C'                      % c...globally standardized
            Vi = (nObs^2)*Vi/(2*a);
        end        
    end
    
    % ---------------------------------------------------------------------
    % Calculate Moran's I:
    % ---------------------------------------------------------------------
    
    mi.mi(i,:) = (mi.My'*Vi*mi.My)/denom;
    
    % --------------------------------------------------------------------- 
    % Calculate eigenvalues and remove k eigenvalues = 0:
    % ---------------------------------------------------------------------
    
    MVM = M*Vi*M;
    MVM = 0.5*(MVM' + MVM);
    evalue = eig(MVM);                     % eigenvalues of MVM
    evalue = flipud(evalue);
    
    min_mi(i,:) = evalue(nObs);             
    max_mi(i,:) = evalue(1);                
    
    for j = 1:nObs                         % sort eigenvalue spectrum 
        if evalue(j) > 0.0000
            idxpos = j;
        else
            break
        end
    end
    
    tau(df,1) = 0;
    tau(1:idxpos) = evalue(1:idxpos);
    tau((idxpos+1):df) = evalue((idxpos+1+k):nObs);
    
    % ---------------------------------------------------------------------
    % Moments of the distribution:
    % --------------------------------------------------------------------- 
    
    % Expected value
    mi.moments_exp(i,:) = sum(tau)/df;   
    
    for j = 1:df
        tauex1(j,:) = (tau(j) - mi.moments_exp(i))^2;
    end
    
    % Variance    
    mi.moments_var(i,:) = (2*sum(tauex1))/(df*(df+2));
    
    
    % Skewness
    for j = 1:df
        tauex1(j,:) = (tau(j) - mi.moments_exp(i))^3;
    end
    mi.moments_skw(i,:) = (8*sum(tauex1))/(df*(df+2)*(df+4));   
    mi.moments_skw(i,:) = mi.moments_skw(i,:)/(mi.moments_var(i,:)^1.5);
    
    % Kurtosis
    for j = 1:df
        tauex1(j,:) = (tau(j) - mi.moments_exp(i))^4;
    end
    for j = 1:df
        tauex2(j,:) = (tau(j) - mi.moments_exp(i))^2;
    end   
    mi.moments_kur(i,:) = (48*sum(tauex1) + 12*(sum(tauex2)^2))/...
        (df*(df+2)*(df+4)*(df+6));
    mi.moments_kur(i,:) = mi.moments_kur(i,:)/(mi.moments_var(i,:)^2);
    
    % ---------------------------------------------------------------------
    % Probability with normal distribution:
    % --------------------------------------------------------------------- 
    
    % z-standardized MI
    mi.z_mi(i,:) = (mi.mi(i)-mi.moments_exp(i))/(sqrt(mi.moments_var(i)));
    % probability
    mi.prob_nv(i,:) = normcdf((mi.mi(i) - mi.moments_exp(i))/...
        (sqrt(mi.moments_var(i))));
    
    if Sad == 'j'
        
        % -----------------------------------------------------------------
        % Saddlepoint approximation:
        % -----------------------------------------------------------------
        
        for j = 1:df
            taumi(j,:) = tau(j) - mi.mi(i);
        end
        
        if i == 1
            
            % -------------------------------------------------------------
            % Secant-search method for global MI:
            % -------------------------------------------------------------
            
            l = 1/(2*taumi(df)) + 0.01;
            h = 1/(2*taumi(1)) - 0.01;
            
            j = 1:df;
            flv(j) = taumi(j)./(1-2*l-taumi(j));
            fl = sum(flv);
            fhv(j) = taumi(j)./(1-2*h-taumi(j));
            fh = sum(fhv);
            
            if fl < 0
                xl = l;
                xh = h;  
            else
                xl = h;
                xh = l;
                swap = fl;
                fl = fh;
                fh = swap;
            end
            
            dx = xh - xl;
            del = 1;
            f = 1;
            z = 10000;                    % number of possible iterations
            for q = 1:z                   % (increase if necessary)  
                
                rtf = xl + dx*fl/(fl - fh);
                fv(j) = taumi(j)./(1-2*rtf*taumi(j));
                f = sum(fv);
                
                if f < 0 
                    del = xl - rtf;
                    xl = rtf;
                    fl = f;
                else
                    del = xh - rtf;
                    xh = rtf;
                    fh = f;
                end
                dx = xh - xl;
                fehler = abs(f);
                
                % displays number of iterations and current error
                % (takes a lot of computation time)                               
                % clc;
                % fprintf('Calculation of global saddlepoint:\n');
                % fprintf('Iteration %g of %g, Error: %g\n',q,noglsad,fehler);
                if abs(f) < 0.00001;        % convergence criterion
                    break;
                end
            end
            mi.sad(i,:) = rtf;
            
            
        else
            
            % -------------------------------------------------------------
            % Exact saddlepoint (root) for local MI:
            % -------------------------------------------------------------
            
            l = tau(df);
            h = tau(1);
            n = df-2;
            aroot = n*mi.mi(i)*(l + h - 2*mi.mi(i)) + mi.mi(i)*(3*l + ...
                    3*h - 4*mi.mi(i)) - 2*l*h;
            broot = (n + 2)*mi.mi(i)*(l - mi.mi(i))*(h - mi.mi(i));
            c1root = l^2 * mi.mi(i)^2 * (n + 1)^2 + h^2 * mi.mi(i)^2*...
                     (n + 1)^2;
            c2root = 2*l*h * (2*l*h - 2*l*mi.mi(i) - 2*h*mi.mi(i) ...
                     - 2*n*mi.mi(i)^2 - n^2 * mi.mi(i)^2 + mi.mi(i)^2);
            mi.sad(i,:) =  0.25*((aroot - sqrt(c1root + c2root))/broot);
            
        end
        
        % -----------------------------------------------------------------
        % Calculation of r and u:
        % -----------------------------------------------------------------
        
        glsad = isnan(mi.sad(i));
        switch glsad            
            
            case 1              % no global saddlepoint was found
                mi.sad_r(i,:) = NaN;
                mi.sad_u(i,:) = NaN;
                mi.prob_sad(i,:) = NaN;
            otherwise           % global saddlepoint was found
                if mi.sad(i) < 0 
                    for j = 1:df
                        r1(j,:) = log(1-2*mi.sad(i)*taumi(j));
                    end
                    mi.sad_r(i,:) = -sqrt(sum(r1));
                end
                if mi.sad(i) > 0 
                    for j = 1:df
                        r2(j,:) = log(1-2*mi.sad(i)*taumi(j));
                    end
                    mi.sad_r(i,:) = sqrt(sum(r2));  
                end
                for j = 1:df
                    uv(j,:) = (taumi(j)^2)/(1-2*mi.sad(i)*taumi(j))^2;
                end
                mi.sad_u(i,:) = mi.sad(i)*sqrt(2*sum(uv));
                
                probsadok = isreal(mi.sad_r(i)) + isreal(mi.sad_u(i));
                if probsadok == 2
                    mi.prob_sad(i,:) = normcdf(mi.sad_r(i)-(1/mi.sad_r(i))...
                                       *log(mi.sad_r(i)/mi.sad_u(i)),0,1);
                else
                    mi.sad_r(i,:) = NaN;
                    mi.sad_u(i,:) = NaN;
                    mi.prob_sad(i,:) = NaN;
                end                   
        end
    end
end

t = toc;        % stop timer
clc

% -------------------------------------------------------------------------
% Display results:
% -------------------------------------------------------------------------

mi.rechenzeit = ['The calculations took ',num2str(t),' seconds.'];

if Sad == 'j'
    if q == z
        mi.glsadinfo = ['It was not possible to find a global saddlepoint.'];
    else
        mi.glsadinfo = ['The global saddlepoint was found after ',num2str(q) ...
                        ' iterations.'];
    end
end

