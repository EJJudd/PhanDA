function d18OpH = pHCorrect(d18oc,T,S,pH,type)
% function d18OpH = pHCorrect(d18oc,T,S,pH,type)
%
% this is an add-on for the BAYFOX forward model that corrects
% forward-modeled d18Oc for the "pH effect", sometimes called the
% "carbonate ion effect". The BAYFOX calibration uses coretop sediments
% that were deposited at presumably preindustrial/late Holocene ocean pH.
% However ocean pH changes with atmospheric pCO2 so for deep-time
% applications, a correction needs to be made to account for this. The pH
% of the ocean influences the speciation of dissolved inorganic carbon,
% thus impacting the net fractionation factor when calcite is formed.
%
% ----- Inputs -----
% d18Oc: The 1000-member ensemble estimate of d18O of calcite from bayfox_forward (N x 1000)
%
% T: SST (scalar or vector)
% S: SSS (scalar or vector)
% pH: pH values for the paleo time period (scalar)
% type: Optional. Enter 1 to use the reduced sensitivity of O. universa.
%
% ----- Outputs -----
%
% d18opH: A 1000-member ensemble estimate of d18O of calcite, corrected for
% the pH effect.
%

% process optional call to a reduced sensitivity
ng=nargin;
if ng == 5
else
    type = 0;
end
%
%ensure vector for pH
pH = pH(:);
%normalize the curve to PI pH
PIpH = 8.16; %best estimate of preindustrial pH. 
%SOURCE: Jiang et al., 2019, https://doi.org/10.1038/s41598-019-55039-4,
%1770-1850 global mean (weighted by latitude).
% set pH range for theoretical curve
% pHrange = (6.5:.01:8.8)';
% use alphaT to get the fractionation factor alpha
AT = alphaT(T,S,pH);
ATb = alphaT(T,S,PIpH);
% calculate carbonate and convert to VPDB scale
% don't need to account for d18Osw since we are interested in relative change.
d18Ocf = ((AT - 1).*1000).*0.97001 - 29.99; %Brand et al 2014 conversion
d18Ocb = ((ATb - 1).*1000).*0.97001 - 29.99;
%normalize data
d18OcN = d18Ocf - d18Ocb;
if type == 1
    %use 1-sigma for Orbulina
    dcErr = 0.14;
    %Orbulina regression
    d18OpH = d18oc + normrnd(repmat(-0.27.*(pH - PIpH),size(d18oc,1),size(d18oc,2)),repmat(dcErr,size(d18oc,1),size(d18oc,2)));
else
    %define the 1-sigma error. Best estimate based on culture studies is 0.27
    dcErr = 0.27;
    %normal random sampling
    d18OpH = d18oc + normrnd(repmat(d18OcN,1,size(d18oc,2)),repmat(dcErr,size(d18oc,1),size(d18oc,2)));
end

end