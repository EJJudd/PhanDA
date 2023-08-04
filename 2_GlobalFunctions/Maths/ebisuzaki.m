function [F,rReal,rSim,cv] = ebisuzaki(x,y,sig,nsim)
 
% EBISUZAKI - a method for estimating significant correlation values for
%             autocorrelated time series
%
%   [F] = ebisuzaki(x,y,sig,nsim)
%
% This function creates 'nsim' random time series that have the same power
% spectrum as the original time series but with random phases.
%
% Input
% x : vector of (real) numbers with the first time series 
% y : vector of (real) numbers with the second time series
% sig : [optional] significance level for critical value estimation [0.05]  
% nsim : [optional] the number of simulations [10000]
%
% Output
% F : Fraction of time series with higher correlation coefficents
%
% Ebisuzaki, W, 1997: A method to estimate the statistical 
% significance of a correlation when the data are serially correlated.  
% J. of Climate, 10, 2147-2153.
%
% modified by Kevin Anchukaitis from the original Matlab function by Vincent Moron
% and the original C code by W. Ebisuzaki

% 2008, Kevin Anchukaitis, Lamont Doherty Earth Observatory, kja@ldeo.columbia.edu
% Ebisuzaki's C code: ftp://ftp.cpc.ncep.noaa.gov/wd51we/random_phase
 
if nargin < 4; nsim = 10000; end
if nargin < 3; sig = 0.05; end
if nargin < 2; error('Minimum input is the 2 time series vectors'); end

rand('seed',-22459); randn('seed',-22459); 
% Exact from Ebisuzaki's C code
% rand('seed',sum(100*clock)); randn('seed',sum(100*clock)); % better?
% rand('twister',sum(100*clock)); randn('twister',sum(100*clock)); % best?

n=length(x);
n2=floor(n/2);
x=x(:); y=y(:);
x=standardize(x);
y=standardize(y);

if size(x,1)~=size(y,1); error('Size of x and y must be the same'); end

[rr,~] = corrcoef(x,y,'rows','pairwise'); rReal = rr(2,1);

xf=fft(x); yf=fft(y);

modx=abs(xf); mody=abs(yf);  
X = NaN(size(x,1),nsim); Y = NaN(size(y,1),nsim);

rSim=NaN(nsim,1);
for i=1:nsim % this is very slow ...
   if n/2==n2
      %an1=angle(xf(2:n2));
      tt=randn(n2-1,1);
      anr1=(tt.*2.*pi)-pi;
      anr1=[0;anr1(:);0;flipud(anr1(:).*(-1))];
      recf=modx.*exp(sqrt(-1)*anr1); 
      X(:,i)=real(ifft(recf));
 
      %an1=angle(yf(2:n2));
      tt=randn(n2-1,1);
      anr1=(tt.*2.*pi)-pi;
      anr1=[0;anr1(:);0;flipud(anr1(:).*(-1))];
      recf=mody.*exp(sqrt(-1)*anr1);
      Y(:,i)=real(ifft(recf));            
      
   else 
      %an1=angle(xf(2:n2+1));
      tt=randn(n2,1);
      anr1=(tt.*2.*pi)-pi;
      anr1=[0;anr1(:);flipud(anr1(:).*(-1))];
      recf=modx.*exp(sqrt(-1)*anr1); 
      X(:,i)=real(ifft(recf));

      %an1=angle(yf(2:n2+1));
      tt=randn(n2,1);
      anr1=(tt.*2.*pi)-pi;
      anr1=[0;anr1(:);flipud(anr1(:).*(-1))];
      recf=mody.*exp(sqrt(-1)*anr1); 
      Y(:,i)=real(ifft(recf));      
      
   end

[rs,~] = corrcoef(standardize(X(:,i)),standardize(Y(:,i))); 
rSim(i) = rs(2,1);

end

F = sum(abs(rSim) > abs(rReal))/nsim;

% Find critical value
rSims = sort(abs(rSim(:))); 
cvs = floor(nsim - (nsim*sig));
cv = rSims(cvs);

end

%% NOTES: Compare Ebisuzaki's code with this one:
 
% bash$ r_phase 0.01 64 x y
% reading x
% reading y
% calculating
% testing x and y (64 points) for 0.010000 significance
% sample |corr| 0.078581,  fraction of samples with larger |corr| 0.519120
% critical |corr| 0.304759 at 0.010000 sig level, 50000 random samples used

% ebisuzaki(xE,yE,0.01,50000);
% Creating synthetic time series...
% Performing Monte Carlo significance test
% ---
% Observed correlation coefficent: 0.078581
% Fraction of |coefficients| larger than Observed: 0.521 [Threshold Value: 0.01]
% Critical R Value: 0.30726


