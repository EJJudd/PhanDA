function d18oc = d18foram_forward(t,d18osw,species,bayes)
% function d18oc = bayfox_forward(t,d18osw,species,bayes)
%
% BAYFOX forward model for d18O of planktic foraminifera.
% Predicts d18O in foraminiferal calcite from given temperature, d18O of
% seawater, and species group.
% ----- Inputs -----
% t: A scalar or vector of sea surface temperatures (Celsius) (1 x N) or (N x 1)
%
% d18osw: A scalar or vector of d18O of seawater (in VSMOW) (1 x N) or (N x 1)
%
% species: A character string corresponding the foram species. Options are:
% 'bulloides' = G. bulloides
% 'incompta' = N. incompta
% 'pachy' = N. pachyderma
% 'ruber' = G. ruber
% 'sacculifer' = T. sacculifer
% 'all' = use the pooled annual (non-species specific) model
% 'all_sea' = use the pooled seasonal (non-species specific) model
%
% bayes: (optional - mainly for use with the DASH interface) A Bayesian
% posterior to use for the calibration. if empty, the default posterior
% file is used.
%
% ----- Outputs -----
%
% d18oc: A 2000-member ensemble estimate of d18O of calcite. (N x 2000)

    % Ensure column vectors.
    t=t(:);
    d18osw=d18osw(:);
    species_list = {'bulloides','incompta','pachy','ruber','sacculifer','all','all_sea'};
    %check that you have an allowable species
    if ~ismember(species,species_list)
        error('Species not recognized');
    end
    % process optional call to a different posterior file
    ng=nargin;
    if ng == 4
    else
        bayes=["poolann_params.mat";"poolsea_params.mat";"hiersea_params.mat"];
    end
    %
    id = (1:1:5);
    %define dimensions
    Nobs=length(t);
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
        id = id(ismember(species_list,species));
    end
    %get posterior params
    betaT=params.betaT(:,id);
    alpha=params.alpha(:,id);
    sigma=params.sigma(:,id);

    % Unit adjustment for permil VSMOW to permil VPDB.
    d18osw_adj = d18osw - 0.27;

    %vectorized calculation of ensemble.
    d18oc = normrnd(alpha' + t * betaT' + d18osw_adj, repmat(sigma',Nobs,1));
end