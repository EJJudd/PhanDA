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
% If working off of older input workspace, manually enter the date
% (Convention: ddMmmyyyy)
filedate = '27Jul2023';
AssDir = [OneDrive,'/AssimilationOutputs/PhanerozoicDA_',filedate];
InputDir = [AssDir,'/InputWorkspaces'];
OutputDir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Code/DataAssimilation/5_NonGlobalFiles/Rexperiments';
FigDir = '/Users/emilyjudd/Library/CloudStorage/OneDrive-SyracuseUniversity/PhanTASTIC/Figures/Supplemental';
savefig = true;
% (c) Load the ensemble, YYe values, exp info, and GTS file
EnsDir = [iCloud,'/ModelOutputs/AssimilationFiles_beta'];
enfilename = [EnsDir,'/HadCM3_all.ens'];
Ens = ensemble(enfilename);
EnsMeta = Ens.metadata;
load([InputDir,'/YYe.mat'])
load("ExpNo.mat", "ExpNo")
load("GTS2020_PETM.mat");

%% PART 1: DECISION MAKING
% (a) Select stage numbers to assimilate
DArange = [1,4,8];
stagelab = strrep( strcat( ...
        'S', num2str(DArange','%02.f'), '_', GTS.Stage(DArange)), ' ', '');
% (b) Define ensemble
% How should the ensemble members be indexed?
dims = string(EnsMeta.ensembleDimensions{1});
ExpNoAll = EnsMeta.ensemble{1,1}{1, dims == "expno"};

% (c) Define R scenarios
Rvals.d18cforam = [.00001 .000025 .00005 .000075 .0001 .00025 .0005 .00075 .001 .0025 .005 .0075 .01 .025 .05 .075 .1 .25 .5 .75 1];
Rvals.tex = [.0000001 .00000025 .0000005 .00000075 .000001 .0000025 .000005 .0000075 .00001 .000025 .00005 .000075 .0001 .00025 .0005 .00075 .001 .0025 .005 .0075 .01 .025 .05 .075 .1];
Rvals.uk = [.0000001 .00000025 .0000005 .00000075 .000001 .0000025 .000005 .0000075 .00001 .000025 .00005 .000075 .0001 .00025 .0005 .00075 .001 .0025 .005 .0075 .01 .025 .05 .075 .1];
Rvals.mg = [.0000025 .000005 .0000075 .00001 .000025 .00005 .000075 .0001 .00025 .0005 .00075 .001 .0025 .005 .0075 .01 .025 .05 .075 .1 .2 .3 .4 .5 .6 .7 .8 .9 1 1.5 2.5];

% Define constants
dimorder = ["lat";"lon"];
AssumptionFiles.predictsw = load('SeawaterMdl.mat');

% Add pH correction for the deep time stages
Ye.S04_Calabrian.d18cforam = Ye.S04_Calabrian.d18cforam + ...
    Assumptions.S04_Calabrian.local.d18cforam.d18OpHrec;
Ye.S06_Piacenzian.d18cforam = Ye.S06_Piacenzian.d18cforam + ...
    Assumptions.S06_Piacenzian.local.d18cforam.d18OpHrec;
Ye.S08_Messinian.d18cforam = Ye.S08_Messinian.d18cforam + ...
    Assumptions.S08_Messinian.local.d18cforam.d18OpHrec;

proxies = ["d18cforam";"mg";"tex";"uk"];
for a = 1:numel(DArange)
    % Parse ensemble to use
    sidx = find(ExpNo.Stage == GTS.Stage(DArange(a)));
    expidx = [ExpNo.ExpNo1(sidx), ExpNo.ExpNo2(sidx)];
    ensidx = any(ExpNoAll == expidx,2);
    ensuse = Ens.useMembers(ensidx);
    ensusemeta = EnsMeta.removeMembers(~ensidx);
    ensusemat = load(ensuse);
    for ii = 1:numel(proxies)
        N = numel(Y.(stagelab{a}).(proxies{ii}));
        if N >= 5
        LOO.(proxies{ii}).(stagelab{a}) = NaN(numel(Rvals.(proxies{ii})),N);
        GMST.(proxies{ii}).(stagelab{a}) = NaN(numel(Rvals.(proxies{ii})),N);
        CE.(proxies{ii}).(stagelab{a}) = NaN(numel(Rvals.(proxies{ii})),1);
        MAE.(proxies{ii}).(stagelab{a}) = NaN(numel(Rvals.(proxies{ii})),1);   
        d18swglobal = Assumptions.(stagelab{a}).global.d18Osw.raw+Assumptions.(stagelab{a}).global.d18Osw.Snowball;
        pHglobal = Assumptions.(stagelab{a}).global.pHrec;
        for jj = 1:numel(Rvals.(proxies{ii}))
            for kk = 1:N
                assimilate = true(N,1);
                assimilate(kk) = false;
                % Run assimilation
                Ruse = repmat(Rvals.(proxies{ii})(jj),N-1,1);
                kf = kalmanFilter;
                kf = kf.prior(ensusemat);
                kf = kf.observations(Y.(stagelab{a}).(proxies{ii})(assimilate));
                kf = kf.uncertainties(Ruse);
                kf = kf.estimates(Ye.(stagelab{a}).(proxies{ii})(assimilate,:));
                output = kf.run;
                % Process results
                tasmean = kel2cel(ensusemeta.regrid("tas", output.Amean, 'order', dimorder));
                GMST.(proxies{ii}).(stagelab{a})(jj,kk) = latweightgmst(tasmean);
                if proxies{ii} == "uk"
                    LOO.(proxies{ii}).(stagelab{a})(jj,kk) = runukPSM_DAresults(...
                        ensuse, ensusemeta, UPD.(stagelab{a}), output, ~assimilate);
                elseif proxies{ii} == "tex"
                    LOO.(proxies{ii}).(stagelab{a})(jj,kk) = runtexPSM_DAresults(...
                        ensuse, ensusemeta, UPD.(stagelab{a}), output, ~assimilate);
                elseif proxies{ii} == "mg"
                     [Stos,Somega,Sso,SpH,Staxon] = runmgcaPSM_DAresults(ensuse, ...
                         ensusemeta, UPD.(stagelab{a}), output,  ~assimilate);
                    % Do calculation outisde of function, because for some
                    % reason it really slows down the script?
                    lnmg = mgcaPSM_forward(UPD.(stagelab{a}).mg.Age(kk), Stos, Somega, Sso, ...
                            SpH, UPD.(stagelab{a}).mg.CleaningMethod(kk), Staxon, 1);
                    LOO.(proxies{ii}).(stagelab{a})(jj,kk) = median(lnmg,2);
                elseif proxies{ii} == "d18cforam"
                    LOO.(proxies{ii}).(stagelab{a})(jj,kk) = rund18cforamPSM_DAresults(...
                        ensuse, ensusemeta, ensusemat, UPD.(stagelab{a}), output, ...
                           d18swglobal, pHglobal, AssumptionFiles, ~assimilate);
                elseif proxies{ii} == "d18acmacro"
                    LOO.(proxies{ii}).(stagelab{a})(jj,kk) = rund18acmacroPSM_DAresults(...
                        ensuse, ensusemeta, ensusemat, UPD.(stagelab{a}), ...
                        output, d18swglobal, pHglobal, AssumptionFiles, ~assimilate);
                elseif proxies{ii} == "d18p"
                    LOO.(proxies{ii}).(stagelab{a})(jj,kk) = rund18pPSM_DAresults(...
                        ensuse, ensusemeta, ensusemat, UPD.(stagelab{a}), ...
                        output, d18swglobal, AssumptionFiles, ~assimilate);
                end
                disp(kk)
            end
            CE.(proxies{ii}).(stagelab{a})(jj,1) = coefefficiency(Y.(stagelab{a}).(proxies{ii}),LOO.(proxies{ii}).(stagelab{a})(jj,:)');
            MAE.(proxies{ii}).(stagelab{a})(jj,1) = sum(abs(Y.(stagelab{a}).(proxies{ii})-LOO.(proxies{ii}).(stagelab{a})(jj,:)'))/N;
        fprintf('Completed Rval %.0f/%.0f (%s)\n', jj, numel(Rvals.(proxies{ii})), proxies{ii})
        end
        end
    end
end
save([OutputDir,'/Rexperiment.mat'],'Rvals','GMST','CE','MAE','LOO')

%% PART 3: PLOT RESULTS
load([OutputDir,'/Rexperiment.mat'])
% Get stage names
DArange = [1,4,6,8];
load("GTS2020_PETM.mat")
stagelab = GTS.Stage(DArange);
agelab = ["(11.7 - 0 ka)";"(1.8-0.774 Ma)";"(3.6-2.58 Ma)";"(7.25-5.33 Ma)"];
% Remove mg/ca Messinian data (negative CE)
CE.mg = rmfield(CE.mg,'S08_Messinian');
MAE.mg = rmfield(MAE.mg,'S08_Messinian');
% Plot
fig = figure('Position',[50 200 1100 450],'Color','w');
proxies = ["d18cforam","mg","tex","uk"];
ax = gobjects(numel(proxies)*2);
t = tiledlayout(2,numel(proxies),'Padding','none','TileSpacing','none');
cm = flipud(hex2rgb({'#9B2226','#CA6702','#FFB703','#0A9396','#005f73'},1));
clear Ruse
Ruse.d18cforam = [1e-3,1e-2,1e-1];
Ruse.tex = [1e-4,1e-3,1e-2];
Ruse.uk = [2.5e-5,2.5e-4,2.5e-3];
Ruse.mg = [5e-3,5e-2,5e-1];
titles = ["δ^{18}O_{carbonate} (forams)","Mg/Ca","TEX_{86}","U^{K'}_{37}"];

for ii = 1:numel(proxies)
    ax(ii) = nexttile(ii); hold on, box on
    ax(ii+numel(proxies)) = nexttile(ii+numel(proxies)); hold on, box on
    stages = sort(fieldnames(CE.(proxies{ii})));
    rectangle(ax(ii),'Position',[Ruse.(proxies{ii})(1),0,...
        range(Ruse.(proxies{ii})),5],'FaceColor',[.5 .5 .5 .5],'EdgeColor','none')
    rectangle(ax(ii+numel(proxies)),'Position',[Ruse.(proxies{ii})(1),-4,...
        range(Ruse.(proxies{ii})),5],'FaceColor',[.5 .5 .5 .5],'EdgeColor','none')
    for jj = 1:numel(Ruse.(proxies{ii}))
        plot(ax(ii),[Ruse.(proxies{ii})(jj),Ruse.(proxies{ii})(jj)],...
            [0 5],'-','LineWidth',1.5,'Color',[.5 .5 .5])
        plot(ax(ii+numel(proxies)),[Ruse.(proxies{ii})(jj),Ruse.(proxies{ii})(jj)],...
            [-4 1],'-','LineWidth',1.5,'Color',[.5 .5 .5])
    end
    for jj = 1:numel(stages)
        plot(ax(ii),Rvals.(proxies{ii}),MAE.(proxies{ii}).(stages{jj}),...
            'o-','Color',cm(jj,:),'MarkerFaceColor',cm(jj,:))
        plot(ax(ii+numel(proxies)),Rvals.(proxies{ii}),CE.(proxies{ii}).(stages{jj}),...
            'o-','Color',cm(jj,:),'MarkerFaceColor',cm(jj,:))
    end
    set(ax(ii),'XScale',"log",'XTick',[10^-8,10^-7,10^-6,10^-5,10^-4,10^-3,...
        10^-2,10^-1,10^0],'XLim',[min(Rvals.(proxies{ii})),max(Rvals.(proxies{ii}))])
    set(ax(ii+numel(proxies)),'XScale',"log",'XTick',[10^-8,10^-7,10^-6,10^-5,...
        10^-4,10^-3,10^-2,10^-1,10^0],'XLim',[min(Rvals.(proxies{ii})),max(Rvals.(proxies{ii}))])
    title(ax(ii),titles(ii),'FontWeight','bold','FontName','Arial','FontSize',13,'Color','k')
end
% Manually adjust ylimits
ylim(ax(1),[.15 1])
ylim(ax(2),[.1 .25])
ylim(ax(3),[0.04 .08])
ylim(ax(4),[.04 .15])
ylim(ax(5),[0 1])
ylim(ax(6),[0 .95])
ylim(ax(7),[0 1])
ylim(ax(8),[.4 1])
xlabel(t,'R value','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(ax(1),'Mean absolute error','FontName','Arial','FontSize',13,'FontWeight','bold')
ylabel(ax(numel(proxies)+1),'Coefficient of efficiency','FontName','Arial','FontSize',13,'FontWeight','bold')    

text(ax(5),1.5e-5,.55,stagelab(1),'FontName','Arial','FontSize',13,'FontWeight','bold','Color',cm(1,:))
text(ax(5),1.5e-5,.485,agelab(1),'FontName','Arial','FontSize',11,'FontWeight','normal','Color',cm(1,:))

text(ax(5),1.5e-5,.41,stagelab(2),'FontName','Arial','FontSize',13,'FontWeight','bold','Color',cm(2,:))
text(ax(5),1.5e-5,.345,agelab(2),'FontName','Arial','FontSize',11,'FontWeight','normal','Color',cm(2,:))

text(ax(5),1.5e-5,.27,stagelab(3),'FontName','Arial','FontSize',13,'FontWeight','bold','Color',cm(3,:))
text(ax(5),1.5e-5,.205,agelab(3),'FontName','Arial','FontSize',11,'FontWeight','normal','Color',cm(3,:))

text(ax(5),1.5e-5,.13,stagelab(4),'FontName','Arial','FontSize',13,'FontWeight','bold','Color',cm(4,:))
text(ax(5),1.5e-5,.065,agelab(4),'FontName','Arial','FontSize',11,'FontWeight','normal','Color',cm(4,:))


if savefig
    export_fig(gcf,[FigDir,'/SupFig_Rexperiments'],'-p0.01','-m5')
end
