function [lmg, mgsw] = ...
    mgcaPSM_forward(age,t,omega,salinity,pH,clean,species,varargin)
% function lmg=baymag_forward_ln(age,t,omega,salinity,pH,clean,species,varargin)
%
% BAYMAG forward model for Mg/Ca of planktic foraminifera.
% Predicts Mg/Ca in foraminiferal calcite from given temperature, salinity,
% omega, cleaning method, and species group.
%
% ----- Input -----
% age = scalar or N x 1 vector of ages. Needs to be Ma for seawater
% correction! otherwise units don't matter. 
% t = scalar or N x 1 vector of SST values.
% omega = scalar or N x 1 vector of bottom water saturation state
% salinity = scalar or N x 1 vector of salinity (psu)
% pH = scalar or N x 1 vector of pH (total scale). If you are using a
% species not sensitive to pH you can enter a dummy value.
% clean = scalar to describe cleaning technique. options:
%   1 = reductive 
%   0 = oxidative
%   values between 0 and 1 are allowed and will be treated as a mix of
%   cleaning methods.
% species = string of target species. six options:
%   'ruber' = Globigerinoides ruber, white or pink
%   'bulloides' = Globigerina bulloides
%   'sacculifer' = Trilobatus sacculifer
%   'pachy' = Neogloboquadrina pachyderma
%   'incompta' = Neogloboquadrina incompta
%   'all' = pooled calibration, annual SST
%   'all_sea' = pooled calibration, seasonal SST
% varargin = three optional arguments:
%  1: a scalar to choose whether to account for changes in mgca of seawater.
%  For this to work properly your ages need to be in units of *millions of years*
%  (Ma)! If this is not entered then no seawater correction is applied.
%    0 = do not include a sw term
%    1 = include a sw term
%  2 and 3: a prior mean and prior standard deviation, scalar values, in
%  Mg/Ca units of mmol/mol. If not entered, it is assumed that no prior is
%  required.
%  4: input different bayesian parameters. Should be a string array of size
%  3 x 1, with each column containing the parameters filename.
%
% ----- Output -----
% mg = N x 1000 ensemble of Mg/Ca values
%
% ----- Dependencies -----
% 1) pooled_model_params.mat, pooled_sea_model_params.mat,
% species_model_params.mat: .mat files with default Bayesian parameters
% ----- Examples -----
%
% mg = baymag_forward(age,t,1.05,35,8,0,'ruber');
% mg = baymag_forward(age,t,1.05,35,8,0,'ruber',1);
% mg = baymag_forward(age,t,1.05,35,8,0,'ruber',1,4,2);
%
% function created by Dr. Jessica Tierney, The University of Arizona (2019)

%% deal with optional arguments
ng=nargin;
if ng==11
    sw=varargin{1};
    prior_mean=varargin{2};
    prior_sig=varargin{3};
    bayes=varargin{4};
elseif ng==10
    sw=varargin{1};
    prior_mean=varargin{2};
    prior_sig=varargin{3};
    bayes=["pooled_model_params.mat";"pooled_sea_model_params.mat";"species_model_params.mat"];
elseif ng==9
    error('You entered a prior mean, but not a prior standard deviation');
elseif ng==8
    sw=varargin{1};
    bayes=["pooled_model_params.mat";"pooled_sea_model_params.mat";"species_model_params.mat"];
elseif ng==7
    sw=0;
    bayes=["pooled_model_params.mat";"pooled_sea_model_params.mat";"species_model_params.mat"];
else
    error('You entered too many or too few arguments');
end
%% prepare data and parameters
%ensure everything is column vectors.
t=t(:);
salinity=salinity(:);
omega=omega(:);
pH=pH(:);
clean=clean(:);
species_list = {'ruber','bulloides','sacculifer','pachy','incompta','all','all_sea'};
species_list_model = {'ruber','bulloides','sacculifer','pachy'};
%check that you have an allowable species
if ~ismember(species,species_list)
    error('Species not recognized');
end
%check for incompta. Change to pachy.
if ismember(species,{'incompta'})
    species='pachy';
else
end 
id = (1:1:4);
%define dimensions
Nobs=length(t);
%assign variables and vectorize
omega=(omega.^-2).*ones(Nobs,1);
clean=clean.*ones(Nobs,1);
salinity=salinity.*ones(Nobs,1);
pH=pH.*ones(Nobs,1);
%load appropriate model
if  strcmp(species,'all')
    params = load(bayes(1));
    id = 1;
elseif strcmp(species,'all_sea')
    params = load(bayes(2));
    id = 1;
else
    params = load(bayes(3));
    %grab id location for species.
    id = id(ismember(species_list_model,species));
end
%get posterior params
betaT=params.betaT;
betaC=params.betaC;
betaO=params.betaO;
betaS=params.betaS;
betaP=params.betaP;
%for sigma and alpha, grab the appropriate species
sigma=params.sigma(:,id);
alpha=params.alpha(:,id);

Nparams=length(betaT);

%get mgca sw estimates
if sw==1
   %load MgCa_sw parameters
    load('mgsw_iters.mat','xt','mg_smooth');
    %hold onto modern value
    mg_mod=mg_smooth(1,:);
    %interpolate to given ages
    mgsw=interp1(xt,mg_smooth,age);
    %ratio to modern value, take log value
    mgsw=log(mgsw./repmat(mg_mod,Nobs,1));
else
    mgsw=0;
end

%convert prior mean and sig
if ng==10
    prior_mean_log=log(prior_mean);
    prior_sig_log=mean([log(prior_mean + prior_sig)-log(prior_mean) log(prior_mean)-log(prior_mean - prior_sig)]);
else
end

%calculate mg for all params
if id < 3 %if you selected ruber, bulloides, or the "all" models, assume a pH sensitivity.
    lmg_mean=repmat(alpha',Nobs,1) + t * betaT' + omega * betaO' + salinity * betaS' + pH * betaP' + (1-clean * betaC') + mgsw;
    lmg_sig=repmat(sigma',Nobs,1);
    if ng==10 %if you are using a prior, then calculate posterior values.
        post_mean_num = repmat(prior_sig_log,Nobs,Nparams).^-2 .* repmat(prior_mean_log,Nobs,Nparams) + lmg_sig.^-2 .* lmg_mean;
        post_mean_den = repmat(prior_sig_log,Nobs,Nparams).^-2 + lmg_sig.^-2;
        post_mean = post_mean_num ./ post_mean_den;
        post_sig = sqrt(post_mean_den.^-1);
        lmg=post_mean + randn(Nobs,Nparams).*post_sig;
    else
        lmg=lmg_mean + randn(Nobs,Nparams).*lmg_sig;
    end
else %if you selected pachyderma or sacculifer, assume no pH sensitivity
    lmg_mean=repmat(alpha',Nobs,1) + t * betaT' + omega * betaO' + salinity * betaS'+ (1-clean * betaC') + mgsw;
    lmg_sig=repmat(sigma',Nobs,1);
    if ng==10 %if you are using a prior, then calculate posterior values.
        post_mean_num = repmat(prior_sig_log,Nobs,Nparams).^-2 .* repmat(prior_mean_log,Nobs,Nparams) + lmg_sig.^-2 .* lmg_mean;
        post_mean_den = repmat(prior_sig_log,Nobs,Nparams).^-2 + lmg_sig.^-2;
        post_mean = post_mean_num ./ post_mean_den;
        post_sig = sqrt(post_mean_den.^-1);
        lmg=post_mean + randn(Nobs,Nparams).*post_sig;
    else
        lmg=lmg_mean + randn(Nobs,Nparams).*lmg_sig;
    end
end