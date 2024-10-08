%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ASSIMILATION STEP 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Run the Data Assimilation (!!)  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Created: 05/21 (E. Judd)
% Last updated: 08/23 (E. Judd)

% Notes: 

%   -->  This script is compatable with the most recent Dash release
%        (https://github.com/JonKing93/DASH/releases/tag/v4.2.2)


%% PART 1: DEFINE THE ASSIMILATION DIRECTORIES

% ----- UPDATE TO MATCH YOU FILEPATHS AND FOLDER CONFIGURATION -----
% Root directories (in my case, I store the input data -- i.e., the proxy 
% data and model priors -- on iCloud, and save the assimilation outputs on
% OneDrive)
iCloud = '/Users/emilyjudd/Library/Mobile Documents/com~apple~CloudDocs';
OutputDrive = '/Users/emilyjudd/Documents/PhanDA/5_Outputs/';

% Ensemble Directory: (stored on iCloud)
EnsDir = [iCloud,'/ModelOutputs/AssimilationFiles_beta'];
enfilename = [EnsDir,'/HadCM3_all.ens'];
Ens = ensemble(enfilename);
EnsMeta = Ens.metadata;

% Assimilaion Directory: (soted on OneDrive)
ad = '21May2024';
AssDir = [OutputDrive,'/AssimilationOutputs/PhanerozoicDA_',ad];
InputDir = [AssDir,'/InputWorkspaces'];
OutputDir = [AssDir,'/OutputWorkspaces'];
% FigDir = [AssDir,'/OutputFigs'];
% mkdir([FigDir,'/EnsSummary'])
% mkdir([FigDir,'/StageSummary'])

% YYe data:
load([InputDir,'/YYe.mat'])

% ----- DO NOT UPDATE PAST THIS LINE -----

% Other global files:
load("ExpNo.mat", "ExpNo")
load("GTS2020_PETM.mat");

% Define variables necessary for calculations
load("HadCM3Coordinates.mat","Lat")


%% PART 2: RUN DA
% (a) Select stage numbers to assimilate
%     (Cenozoic = 1:23, Mesozoic = 24:53, Paleozoic = 54:101)
DArange = 1:91; 
stagelab = strrep( strcat( ...
        'S', num2str(DArange','%02.f'), '_', GTS.Stage(DArange)), ' ', '');
    
% (b) Define R values    
Rvals.d18cforam = [1e-3,1e-2,1e-1];
Rvals.d18acmacro = [5e-2,.275,5e-1];
Rvals.d18p = [5e-2,.275,5e-1; ...
              .1,.5,.9];
Rvals.tex = [1e-4,1e-3,1e-2];
Rvals.uk = [2.5e-5,2.5e-4,2.5e-3];
Rvals.mg = [5e-3,5e-2,5e-1; ...
                .1 .5 .9];

% (c) Define assimilation options
% Define Rmethods
Ropts = ["low","medium","high"];
% Ropts = "percentile";
% Deine seawater d18O correction
SWopts = ["snowball","veizer","off"];
% SBopts = true;
% Define pH correction
pHopts = ["ens","rec","off"];
% pHopts = true;
% Want to do covariance analysis?
covanalysis = false;

% (d) Define ensemble
% How should the ensemble members be indexed?
dims = string(EnsMeta.ensembleDimensions{1});
ExpNoAll = EnsMeta.ensemble{1,1}{1, dims == "expno"};

% Which variables should be assimilate?
assvar = ["tas","tos"];
allvar = EnsMeta.variables;
noassvar = setdiff(allvar,assvar);

% (e) Preallocate
Ndata = NaN(DArange(end), numel(Ropts)*numel(SWopts)*numel(pHopts));
Nandata = NaN(DArange(end), numel(Ropts)*numel(SWopts)*numel(pHopts));
GMST = cell(DArange(end), numel(Ropts)*numel(SWopts)*numel(pHopts));
LTG = cell(DArange(end), numel(Ropts)*numel(SWopts)*numel(pHopts));
LTGsst = cell(DArange(end), numel(Ropts)*numel(SWopts)*numel(pHopts));
TASpost = cell(DArange(end), numel(Ropts)*numel(SWopts)*numel(pHopts));
TASprior = cell(DArange(end), 1);
ItName = string(zeros(numel(Ropts)*numel(SWopts)*numel(pHopts),1));

c = 1;
% Assimilate through stages
for a = DArange
    % (a) Parse ensemble to use
    sidx = find(ExpNo.Stage == GTS.Stage(a));
    expidx = [ExpNo.ExpNo1(sidx), ExpNo.ExpNo2(sidx)];
    ensidx = any(ExpNoAll == expidx,2);
    ensuse = Ens.useMembers(ensidx);
    ensuse = ensuse.useVariables(assvar);
    ensusemat = load(ensuse);
    ensusemeta = EnsMeta.removeMembers(~ensidx);
    ensusemeta = ensusemeta.remove(noassvar);
    exprunall = ensusemeta.ensemble{1,1}{1, dims == "exprun"};

    if c == 1
        TASprior{a,1} = kel2cel(ensusemeta.regrid("tas",ensusemat,'order',["lat";"lon"]));
    end
    
    % (b) Select correction preferences
    % pH Correction
    for pH = 1:numel(pHopts)
        pHcorr = pHopts(pH);
    % seawater d18O correction
    for SW = 1:numel(SWopts)
        seawatercorr = SWopts(SW);      
    % Rmethod
    for Rmeth = 1:numel(Ropts)
        Rmethod = Ropts(Rmeth);
    
    ItName(c) = sprintf("phCorr = %s\nSeawaterCorr = %s\nRmethod = %s", ...
                     pHcorr, seawatercorr, Rmethod);

    % (c) Assemble the Y, Ye, & R values
    [Yuse, Yeuse, Ruse, proxytype, paleolat, paleolon] = assembleYYeR( UPD.(stagelab{a}), ...
        Y.(stagelab{a}), Ye.(stagelab{a}), Rvals, Assumptions.(stagelab{a}), ...
        pHcorr, seawatercorr, Rmethod, a);
    Ndata(a,c) = numel(Yuse);

    % (c.2) Remove NaNs, but document how many
    %       (Remove later, but use for now)
    nanidx = any(isnan(Yeuse),2);
    Yuse(nanidx) = [];
    Yeuse(nanidx,:) = [];
    Ruse(nanidx) = [];
    Nandata(a,c) = sum(nanidx);
    proxytype(nanidx) = [];
    % (c) Create prior resampling indices
    % [idx, mdls] = randsampensemble(ensusemeta, expidx, newidx);
    [Index.(stagelab{a}){c,1}, exprun, exps, explabs] = indexensemble(ensusemeta,a);

    % (d) Run assimilation
    [gmst, ltg, ltgsst, taspost, ~, Index.(stagelab{a}){c,1}] = ...
        runkf(ensusemeta, ensusemat, Yuse, Ruse, Yeuse, ...
        Index.(stagelab{a}){c,1}, assvar, covanalysis);
    [GMST{a,c},LTG{a,c},LTGsst{a,c},TASpost{a,c}] = rmvthreshold(...
        gmst,ltg,ltgsst,taspost,Index.(stagelab{a}){c,1});
    
    % (f) Convert TASpost to median 
    TASpost{a,c} = median(TASpost{a,c},3);
    

        if Rmeth ~= numel(Ropts)
            c = c+1;
        end
    end
        if SW ~= numel(SWopts)
            c = c+1;
        end
    end
        if pH ~= numel(pHopts)
            c = c+1;
        end
    end
    c = 1;
    message = sprintf('Completed %s (Stage %d/%d)', ...
        GTS.Stage(a), a, numel(DArange));
    printprogress(message, true)

end
save(sprintf('%s/Output.mat', OutputDir), 'GMST', 'LTG',...
    'LTGsst','Index',"Ndata","Nandata","ItName","TASprior","TASpost","Rvals")
    