%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ASSIMILATION STEP 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Run the Data Assimilation (!!)  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Created: 05/21 (E. Judd)
% Last updated: 08/22 (E. Judd)

% Notes: 

%   -->  This script is compatable with the most recent Dash release
%        (https://github.com/JonKing93/DASH/releases/tag/v4.0.0-beta-2)


%% PART 0: DEFINE THE ASSIMILATION DIRECTORIES
% (a) Home Directory (root for all other directories)
% Home Directory (root for all other directories)
OneDrive = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC';
iCloud = '/Users/emilyjudd/Library/Mobile Documents/com~apple~CloudDocs';

% (b) Load the rest of the directories off of the "DirectoryInfo" file stored
% in the InputWorkspaces folder
% If working with assimilation input workspaces made the same day:
% load( sprintf( '%s/AssimilationOutputs/PhanerozoicDA_%s/InputWorkspaces/DirectoryInfo.mat', ...
%                      HomeDir, strrep(date,'-','') ) );
% If working off of older input workspace, manually enter the date
% (Convention: ddMmmyyyy)
filedate = '27Jul2023';
AssDir = [OneDrive,'/AssimilationOutputs/PhanerozoicDA_',filedate];
InputDir = [AssDir,'/InputWorkspaces'];
OutputDir = [AssDir,'/OutputWorkspaces'];
FigDir = [AssDir,'/OutputFigs'];
mkdir([FigDir,'/EnsSummary'])
mkdir([FigDir,'/StageSummary'])

% (c) Load the ensemble, YYe values, exp info, and GTS file
EnsDir = [iCloud,'/ModelOutputs/AssimilationFiles_beta'];
enfilename = [EnsDir,'/HadCM3_all.ens'];
Ens = ensemble(enfilename);
EnsMeta = Ens.metadata;
load([InputDir,'/YYe.mat'])
load("ExpNo.mat", "ExpNo")
load("GTS2020_PETM.mat");

% (d) Define variables necessary for calculations
load("HadCM3Coordinates.mat","Lat")

% (e) Mapping contours
% mapdirectory = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Data/Paleogeographies/Scotese_PaleoDEMS';
% mapfilenames = selectmap(ExpNo.PlateAge, mapdirectory);
% [~,maplat,maplon,dem] = scotesecontourlines(mapfilenames,GTS);

%% PART 1: DECISION MAKING
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
% Deine Snowball correction
SBopts = [true,false];
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
Ndata = NaN(DArange(end), numel(Ropts)*numel(SBopts)*numel(pHopts));
Nandata = NaN(DArange(end), numel(Ropts)*numel(SBopts)*numel(pHopts));
GMST = cell(DArange(end), numel(Ropts)*numel(SBopts)*numel(pHopts));
LTG = cell(DArange(end), numel(Ropts)*numel(SBopts)*numel(pHopts));
LTGsst = cell(DArange(end), numel(Ropts)*numel(SBopts)*numel(pHopts));
TASpost = cell(DArange(end), numel(Ropts)*numel(SBopts)*numel(pHopts));
TASprior = cell(DArange(end), 1);
ItName = string(zeros(numel(Ropts)*numel(SBopts)*numel(pHopts),1));

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
    % Snowball correction
    for SB = 1:numel(SBopts)
        snowballcorr = SBopts(SB);      
    % Rmethod
    for Rmeth = 1:numel(Ropts)
        Rmethod = Ropts(Rmeth);
    
    ItName(c) = sprintf("phCorr = %s\nSnowballCorr = %s\nRmethod = %s\n", ...
                     pHcorr, string(snowballcorr), Rmethod);

    % (c) Assemble the Y, Ye, & R values
    [Yuse, Yeuse, Ruse, proxytype, paleolat, paleolon] = assembleYYeR( UPD.(stagelab{a}), ...
        Y.(stagelab{a}), Ye.(stagelab{a}), Rvals, Assumptions.(stagelab{a}), ...
        pHcorr, snowballcorr, Rmethod, a);
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
    
%     % (e) Print Figures
%     if c == 3
%         plotensresults(gmst, ltg, Index.(stagelab{a}){c,1}, ...
%             stagelab{a}, exprun, exps, explabs, exprunall, ...
%             maplat, maplon, dem, TASprior{a,1}, UPD, ...
%             [FigDir,'/EnsSummary'])
%         plotstagesummary(Yuse, Yeuse, proxytype, paleolat, paleolon,...
%             TASpost{a,c}, GMST{a,c}, stagelab(a), ...
%             maplat, maplon, dem, TASprior{a,1},[FigDir,'/StageSummary'])
%     end

    % (f) Convert TASpost to median 
    TASpost{a,c} = median(TASpost{a,c},3);
    

        if Rmeth ~= numel(Ropts)
            c = c+1;
        end
    end
        if SB ~= numel(SBopts)
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


% [taspriormean, taspriorvar, gmstprior] = tasprior(Ens, EnsMeta, DArange, stagelab, true);
% filename = [OutputDir, '/tasprior.mat'];
% save(filename, "taspriormean", "taspriorvar", "gmstprior")
        
     