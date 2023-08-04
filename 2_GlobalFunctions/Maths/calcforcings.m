%% Function to calculate a solar constant adjusted pCO2 value based on input ages and CO2 values.

function [CO2_corrected, dFtot, dFco2, dFsun, dFsa]  = ...
    calcforcings(age, CO2, albedo, varargin)


%% See equations in:
%   + Foster et al., 2017
%     (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5382278/#b5)
%   + Myhre et al., 1998 
%     (https://doi.org/10.1029/98GL01908)

% For simplicity sake, here I've opted to use the simplified equation to
% calculate the change in radiative CO2 forcing relative to today (dFco2)
% from Myhre et al. (1998) rather than the updated equation of Etminan et
% al. (2016), which includes a more complicated "a" term that accounts for
% changes in nitrous oxide concentrations.

if nargin == 3
    lsmask = [];
elseif nargin == 4
    lsmask = varargin{1};
    land = repmat(0.14,size(CO2,1),1);
elseif nargin == 5
    land = varargin{2};
end
    
% (1) Calculate dFsun (i.e., change in solar forcing relative to today)
    % Present day sun:
    Fs = 1368;
    % Age of Earth
    t0 = 4567;
    % Time since formation of Earth
    t = t0-age;
    % Total solar irradiance at each stage:
    Fts = 1./( 1 + .4 * (1 - t./t0)) * Fs;
    % Calculate energy reaching surface
    Sr = ((1 - albedo) .* Fts) ./ 4;
    % Change in solar forcing relative to today
    dFsun = Sr - Sr(1);
    %dFsun = ( (Fts - Fs) .* (1 - albedo) )./ 4;

% (2) Calculate dFco2
    % Forcing associated with a  CO2 doubling:
    % (Sherwood et al., 2020)
    F2x = 4.0;
    % Initial concentrations (i.e., modern reference value)
    c0 = nanmedian(CO2(1,:));
    % Change in radiative CO2 forcing relative to today
    dFco2 = F2x.*(log(CO2./c0)./log(2));
 
% (3) Calculate land forcing
    % Albedo of land & ocean
    if ~isempty(lsmask)
        % Calculate land area
        L = NaN(size(lsmask,1),1);
        for ii = 1:size(lsmask,1)
            % Land = 1; Ocean = 0;
            L(ii) = latweightgmst(mean(lsmask{ii},3));
        end
        %land = 0.14;
        ocean = 0.06;
        Dac = land - ocean;
        dFsa = -.25*(L-L(1)).*Dac*Fs(1);   

% OLD METHOD:
%         % Latitudinal weight
%         load("HadCM3Coordinates.mat", "Lat")
%         latweight = cosd(Lat);
%         % Latweighted weighted incoming radiation
%         latweightfs = Fs/4*latweight;
%         %latweightfs = latweight;
%         % Calculations:
%         surfalb = lsmask;
%         Fsa = NaN(size(surfalb,1),1);
%         for ii = 1:size(surfalb,1)
%             % calculate the mean latitudinal abedo gradient
%             surfalb{ii}(surfalb{ii} == 1) = land(ii);
%             surfalb{ii}(surfalb{ii} == 0) = ocean;
%             surfalb{ii} = mean(surfalb{ii},[3,2]);
%             % calculate the mean latitudinal forcing gradient
%             surfalbfor = (1-surfalb{ii}).*latweightfs;
%             % calculate the latitudinally weighted surface albedo forcing
%             Fsa(ii) = sum(surfalbfor .* latweight) ./ sum(latweight);
%         end
%         %Express SurfAlb as a change from today
%         dFsa = Fsa-Fsa(1);

    else
        dFsa = 0;
    end
    
% (4) Add the forcings together
    dFtot = dFco2+dFsun+dFsa;

% (5) Convert forcing back into solar adjusted CO2
    CO2_corrected = exp(dFtot*log(2)/F2x)*c0;
end
