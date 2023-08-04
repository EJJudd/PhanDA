function [uk] = ukPSM_forward(ssts, bayes)
% Models Uk'37 from SSTs using the BAYSPLINE calibration.
%
% uk = UK_forward_model( ssts, bayes )
%
% ----- Inputs -----
%
% ssts: A vector of sea surface temperatures (Celsius) (1 x N) or (N x 1)
%
% bayes: (optional - mainly for use with the DASH interface) A Bayesian
% posterior to use for the calibration. if empty, use the default posterior
% file.
%
% ----- Outputs -----
%
% uk: A set of 1000 uk'37 estimates for each SST. (N x 1000)

% Convert ssts to column vector.
ssts=ssts(:);
% process optional call to a different posterior file
% load posteriors for B coefficients and tau^2 variance
ng=nargin;
if ng == 2
    bayesParams=load(bayes);
else
    bayesParams=load('bayes_posterior_v2.mat');
end

%NOTE: calibration is seasonal for the following regions: North Pacific (Jun-Aug),
%North Atlantic (Aug-Oct), Mediterrenean (Nov-May). If the data lie in the
%following polygons then you should provide seasonal SSTs:

%set up spline parameters with set knots
order = 3;%spline order, 3 for quadratic
kn = augknt(bayesParams.knots,order); %knots

%spmak assembles the b-spline with the given knots and current coeffs
bs_b=spmak(kn,bayesParams.bdraws);

%fnxtr linearly extrapolates the spline to evaluate values at SSTs
%outside the calibration range (0-30). w/o this, B-spline will return a NaN
%at SSTs out of this range.
bs=fnxtr(bs_b);
%evaluate the mean value of the spline for your SST obs:
mean_now=fnval(bs,ssts);

%draw from the distribution:
uk=normrnd(mean_now,repmat(sqrt(bayesParams.tau2),1,length(ssts)))';
%any uk values outside 0 to 1 are forced to be in that range.
uk(uk>1)=1;
uk(uk<0)=0;
end