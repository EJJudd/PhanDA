function [x, xm, xs] = standardize(x, scale)
%STANDARDIZE   Centers and normalizes data.
%
%    [XC, XM, XS] = STANDARDIZE(X) centers and normalizes the data in X to
%    XC by subtracting the mean XM of each column and dividing each
%    column by its standard deviation XS. If X contains missing
%    values, indicated by NaNs, the mean and standard deviation of X
%    are computed from the available data.
%
%    [XC, XM, XS] = STANDARDIZE(X, SCALE) centers and normalizes the data
%    in X to zero mean and standard deviation SCALE. The column means
%    are returned as XM and the scale factors as XS = std(X) ./ SCALE.
%
%    Constant columns of X are not scaled.  
%
%    See also CENTER, NANMEAN, MEAN, NANSTD, STD.

  error(nargchk(1,2,nargin))          % check number of input arguments 
  
  if nargin < 2 
    scale = 1;
  end
  
  if ndims(x) > 2,  error('X must be vector or 2-D array.'); end
  
  % if x is a vector, make sure it is a row vector
  if length(x)==prod(size(x))         
    x = x(:);                         
  end 
  [m,n]  = size(x);

  % get mean and standard deviation of x
  if any(any(isnan(x)))               % there are missing values in x
    xm   = nanmean(x);
    xs   = nanstd(x) ./ scale;
  else                                % no missing values
    xm   = mean(x);
    xs   = std(x) ./ scale;
  end

  % test for constant columns
  const  = (abs(xs) < eps);
  nconst = ~const;
  if sum(const) ~= 0
    warning('Constant or nearly constant columns not rescaled.');
    xm   = xm .* nconst + 0*const;
    xs   = xs .* nconst + 1*const;
  end
   
  % remove mean and divide by standard deviation
  x      = (x - repmat(xm, m, 1) ) ./ repmat(xs, m, 1);  
end