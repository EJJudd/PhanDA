%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ASSIMILATION STEP 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% MAKE STATE VECTOR ENSEMBLE %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Created: 05/21 (E. Judd)
% Last updated: 08/22 (E. Judd)

% Notes: 
%   -->  This script is compatable with the most recent Dash release
%        (https://github.com/JonKing93/DASH/releases/tag/v4.0.0-beta-1)
%   -->  This script builds a single state vector/ensemble file that houses
%        all HadCM3 models

clearvars;
close all;
clc;

%% BUILD STATE VECTOR
% (1) Define directories and load grid file
% Home directory
homedir = '/Users/emilyjudd/Library/Mobile Documents/com~apple~CloudDocs/ModelOutputs';
% Gridfile directory
gridfiledir = [homedir,'/','AssimilationFiles_beta'];
% Gridfile name
gridfilename = 'HadCM3_all.grid';
% Load gridfile
HadCM3 = gridfile( [gridfiledir,'/',gridfilename] );

% (2) Name and initialize state vector
sv = stateVector('HadCM3_all');

% (3) Specify & add desired variables
vars = ["tas";"tos";"so";"pr";"tosanom";"soanom";"prlnanom";"tosmm";"dist"];
sv = sv.add( vars, repmat( strcat(gridfiledir,"/",gridfilename), numel(vars) ,1) );

% (4) Edit variables
% Define ensemble dimensions
endims = ["expno"; "exprun"; "expmean"];
sv = sv.design( vars, endims,"ens");
% Specify which grid var to use for in each sv var dimension
meanvars = vars(vars ~= "tosmm");
for ii = 1:numel(meanvars)
    varselect = HadCM3.meta.var == meanvars(ii);
    sv = sv.design(meanvars(ii), "var", "state", varselect);
    if meanvars(ii) == "tos"
        sv = sv.design("tosmm", "var", "state", varselect);
    end
end
% Take the annual means of all but tos_mm
meanvars = vars(vars ~= "tosmm");
sv = sv.mean( meanvars, "time", [], "omitnan");


% (5) Build and save ensemble
enname = strcat(gridfiledir, "/HadCM3_all.ens");
ens = sv.build( 'all', 'file', enname, 'sequential', true, 'overwrite', true);


%% Manipulate ensemble
ens = ensemble(enname);
ensMat = ens.load;
ensMeta = ens.metadata;



