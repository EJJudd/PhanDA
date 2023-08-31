function [Ref, Stage, Age, Temp, PhanDA, Lit, Ridx] = loadcenolitdata(GMST,GTS)

% PART 1: LOAD LITERATURE DATA

% Inglis et al.
% https://doi.org/10.5194/cp-16-1953-2020
% Reporting average & 90% CI
Ref{1,1} = "Inglis et al., 2019";
%Late Paleocene (Thanetian)
Stage{1,1}(1,1) = "Thanetian";
Age{1,1}{1,1} = [56,58];
Temp{1,1}{1,1} = [21.3,26.3,29.1];
%PETM
Stage{1,1}(2,1) = "PETM";
Age{1,1}{2,1} = [55.7,56];
Temp{1,1}{2,1} = [25.9,31.6,35.6];
%EECO (Ypressian)
Stage{1,1}(3,1) = "Ypresian";
Age{1,1}{3,1} = [49.1,53.3];
Temp{1,1}{3,1} = [22.2,27.0,30.7];

% O'Brien et al.
% https://doi.org/10.1073/pnas.2003914117
% Reporting range of data-adjusted GMST
Ref{2,1} = "O'Brien et al., 2020";
%Late Eocene (Pri/Bart)
Stage{2,1}(1,1) = "Bartonian/Priabonian";
Age{2,1}{1,1} = [33.9,40];
Temp{2,1}{1,1} = [22.17,23.45,24.33];
%Early Olig (Rupellian)
Stage{2,1}(2,1) = "Rupelian";
Age{2,1}{2,1} = [33,33.9];
Temp{2,1}{2,1} = [22.57,23.97,25.49];
%Middle Olig (Rupellian)
Stage{2,1}(3,1) = "Rupelian";
Age{2,1}{3,1} = [26.5,33];
Temp{2,1}{3,1} = [21.04,22.38,23.78];
%Late Olig (Chattian)
Stage{2,1}(4,1) = "Chattian";
Age{2,1}{4,1} = [25,23.5];
Temp{2,1}{4,1} = [22.16,23.51,24.68];


% Burls et al.
% https://doi.org/10.1029/2020PA004054
% Reporting +/1 2 std
%   --> Row 1: all data (terrestrial mat & sst)
%   --> Row 2: sst-only
Ref{3,1} = "Burls et al., 2021";
%Late Miocene
Stage{3,1}(1,1) = "Tortonian/Messinian";
Age{3,1}{1,1} = [5.33,11.6];
Temp{3,1}{1,1} = [19.3-0.98*2, 19.3, 19.3+0.98*2; ...
              21.95-0.81*2, 21.95, 21.95+0.81*2];
%Early Middle Miocene
Stage{3,1}(2,1) = "Langhian/Serravallian";
Age{3,1}{2,1} = [11.6,14.5];
Temp{3,1}{2,1} = [21.21-0.56*2, 21.21, 21.21+0.56*2; ...
              24.46-0.81*2, 24.46, 24.46+0.81*2];
%MCO
Stage{3,1}(3,1) = "Burdigalian/Langhian";
Age{3,1}{3,1} = [14.5,16.75];
Temp{3,1}{3,1} = [22.93-1.01*2, 22.93, 22.93+1.01*2; ...
              25.47-1.17*2, 25.47, 25.47+1.17*2];
%Early Middle Miocene
Stage{3,1}(4,1) = "Burdigalian";
Age{3,1}{4,1} = [16.75,20];
Temp{3,1}{4,1} = [21.21-0.56*2, 21.21, 21.21+0.56*2; ...
              24.46-0.81*2, 24.46, 24.46+0.81*2];

% Anagnastou et al. 
% Reporting mean +/- 2 sigma
Ref{4,1} = "Anagnastou et al., 2020";
load("GMST_Anagnastou.mat");
N = 1000;
%Rupelian
span = [1:5];
Stage{4,1}(1,1) = "Rupelian";
Age{4,1}{1,1} = [GMST_Anagnastou.Age(span(1)), GMST_Anagnastou.Age(span(end))];
Dist = normrnd( repmat(GMST_Anagnastou.GMST(span), N ,1), ...
    repmat(GMST_Anagnastou.sigma(span)/2, N, 1) );
Temp{4,1}{1,1} = [mean(Dist)-2*std(Dist), mean(Dist), mean(Dist)+2*std(Dist)];
%Priabonian
span = [6:15];
Stage{4,1}(2,1) = "Priabonian";
Age{4,1}{2,1} = [GMST_Anagnastou.Age(span(1)), GMST_Anagnastou.Age(span(end))];
Dist = normrnd( repmat(GMST_Anagnastou.GMST(span), N ,1), ...
    repmat(GMST_Anagnastou.sigma(span)/2, N, 1) );
Temp{4,1}{2,1} = [mean(Dist)-2*std(Dist), mean(Dist), mean(Dist)+2*std(Dist)];
%Bartonian
span = [16:53];
Stage{4,1}(3,1) = "Bartonian";
Age{4,1}{3,1} = [GMST_Anagnastou.Age(span(1)), GMST_Anagnastou.Age(span(end))];
Dist = normrnd( repmat(GMST_Anagnastou.GMST(span), N ,1), ...
    repmat(GMST_Anagnastou.sigma(span)/2, N, 1) );
Temp{4,1}{3,1} = [mean(Dist)-2*std(Dist), mean(Dist), mean(Dist)+2*std(Dist)];
%Lutetian
span = [54:82];
Stage{4,1}(4,1) = "Lutetian";
Age{4,1}{4,1} = [GMST_Anagnastou.Age(span(1)), GMST_Anagnastou.Age(span(end))];
Dist = normrnd( repmat(GMST_Anagnastou.GMST(span), N ,1), ...
    repmat(GMST_Anagnastou.sigma(span)/2, N, 1) );
Temp{4,1}{4,1} = [mean(Dist)-2*std(Dist), mean(Dist), mean(Dist)+2*std(Dist)];
%Ypresian
span = [83:120];
Stage{4,1}(5,1) = "Ypresian";
Age{4,1}{5,1} = [GMST_Anagnastou.Age(span(1)), GMST_Anagnastou.Age(span(end))];
Dist = normrnd( repmat(GMST_Anagnastou.GMST(span), N ,1), ...
    repmat(GMST_Anagnastou.sigma(span)/2, N, 1) );
Temp{4,1}{5,1} = [mean(Dist)-2*std(Dist), mean(Dist), mean(Dist)+2*std(Dist)];
%PETM
span = [121:149];
Stage{4,1}(6,1) = "PETM";
Age{4,1}{6,1} = [GMST_Anagnastou.Age(span(1)), GMST_Anagnastou.Age(span(end))];
Dist = normrnd( repmat(GMST_Anagnastou.GMST(span), N ,1), ...
    repmat(GMST_Anagnastou.sigma(span)/2, N, 1) );
Temp{4,1}{6,1} = [mean(Dist)-2*std(Dist), mean(Dist), mean(Dist)+2*std(Dist)];


% Tierney et al., 2022
% https://doi.org/10.1073/pnas.2205326119
% Reporting 95% CI
Ref{5,1} = "Tierney et al., 2022";
%PETM
Stage{5,1}(1,1) = "PETM";
Age{5,1}{1,1} = [55.7,56];
Temp{5,1}{1,1} = [33.1,34.1,35.5];
%pre-PETM
Stage{5,1}(2,1) = "Thanetian";
Age{5,1}{2,1} = [56,57];
Temp{5,1}{2,1} = [27.5,28.5,30.1];

% Hansen et al., 2013
% https://doi.org/10.1098/rsta.2012.0294
% Reporting 95% CI
Ref{6,1} = "Hansen et al., 2013";
%Holocene
Stage{6,1}(1,1) = "Holocene";
Age{6,1}{1,1} = [0,0.0117];
Temp{6,1}{1,1} = [14.15,14.15,14.15];

% Ring et al., 2022
% https://doi.org/10.1029/2021PA004364
load("GMST_ring.mat","GMST_Ring")
idx = [5,8,11,14];
% Reporting range of values
Ref{7,1} = "Ring et al., 2022";
%Middle Pliocene
Stage{7,1}(1,1) = "Zanclean";
Age{7,1}{1,1} = [GMST_Ring.MiddlePliocene(1),GMST_Ring.MiddlePliocene(2)];
Temp{7,1}{1,1} = 14.4 + [min(GMST_Ring.MiddlePliocene(idx)),...
    GMST_Ring.MiddlePliocene(idx(end)),max(GMST_Ring.MiddlePliocene(idx))];
%Late Miocene
Stage{7,1}(2,1) = "Tortonian";
Age{7,1}{2,1} = [GMST_Ring.LateMiocene(1),GMST_Ring.LateMiocene(2)];
Temp{7,1}{2,1} = 14.4 + [min(GMST_Ring.LateMiocene(idx)),...
    GMST_Ring.LateMiocene(idx(end)),max(GMST_Ring.LateMiocene(idx))];
%Middle Miocene
Stage{7,1}(3,1) = "Langhian/Burdigalian";
Age{7,1}{3,1} = [GMST_Ring.MiddleMiocene(1),GMST_Ring.MiddleMiocene(2)];
Temp{7,1}{3,1} = 14.4 + [min(GMST_Ring.MiddleMiocene(idx)),...
    GMST_Ring.MiddleMiocene(idx(end)),max(GMST_Ring.MiddleMiocene(idx))];
%Early Miocene
Stage{7,1}(4,1) = "Aquitanian";
Age{7,1}{4,1} = [GMST_Ring.EarlyMiocene(1),GMST_Ring.EarlyMiocene(2)];
Temp{7,1}{4,1} = 14.4 + [min(GMST_Ring.EarlyMiocene(idx)),...
    GMST_Ring.EarlyMiocene(idx(end)),max(GMST_Ring.EarlyMiocene(idx))];
%Early Oligocene
Stage{7,1}(5,1) = "Rupelian";
Age{7,1}{5,1} = [GMST_Ring.EarlyOligocene(1),GMST_Ring.EarlyOligocene(2)];
Temp{7,1}{5,1} = 14.4 + [min(GMST_Ring.EarlyOligocene(idx)),...
    GMST_Ring.EarlyOligocene(idx(end)),max(GMST_Ring.EarlyOligocene(idx))];
%Late Eocene
Stage{7,1}(6,1) = "Priabonian";
Age{7,1}{6,1} = [GMST_Ring.LateEocene(1),GMST_Ring.LateEocene(2)];
Temp{7,1}{6,1} = 14.4 + [min(GMST_Ring.LateEocene(idx)),...
    GMST_Ring.LateEocene(idx(end)),max(GMST_Ring.LateEocene(idx))];
%Middle Eocene
Stage{7,1}(7,1) = "Lutetian";
Age{7,1}{7,1} = [GMST_Ring.MiddleEocene(1),GMST_Ring.MiddleEocene(2)];
Temp{7,1}{7,1} = 14.4 + [min(GMST_Ring.MiddleEocene(idx)),...
    GMST_Ring.MiddleEocene(idx(end)),max(GMST_Ring.MiddleEocene(idx))];
%Early Eocene
Stage{7,1}(8,1) = "Ypresian";
Age{7,1}{8,1} = [GMST_Ring.EarlyEocene(1),GMST_Ring.EarlyEocene(2)];
Temp{7,1}{8,1} = 14.4 + [min(GMST_Ring.EarlyEocene(idx)),...
    GMST_Ring.EarlyEocene(idx(end)),max(GMST_Ring.EarlyEocene(idx))];

% Eichenseer & Jones, 2023
% https://doi.org/10.5194/egusphere-2023-1188
% Reporting 95% CI
Ref{8,1} = "Eichenseer & Jones, 2023";
%EECO
Stage{8,1}(1,1) = "Ypresian";
Age{8,1}{1,1} = []; 
Temp{8,1}{1,1} = [26.7,28.7,30.7];


% PART 2: EXTRACT MATHCING PHANDA VALUE
PhanDA = [];
Lit = [];
Ridx = strings(0,0);
for ii = 1:numel(Ref)
    if ~contains(Ref{ii}, "Anagnastou")
    for jj = 1:numel(Stage{ii})
        s = strsplit(Stage{ii}(jj),'/');
        idx = find(contains(GTS.Stage,s));
        gmst = [];
        for kk = 1:numel(idx)
            gmst = [gmst;GMST{idx(kk)}];
        end
        PhanDA(end+1,:) = prctile(gmst,[5,50,95]);
        Lit(end+1,:) = Temp{ii}{jj}(1,:);
        Ridx(end+1,1) = Ref{ii};
        if contains(Ref{ii}, "Burls")
            PhanDA(end+1,:) = prctile(gmst,[5,50,95]);
            Lit(end+1,:) = Temp{ii}{jj}(2,:);
            Ridx(end,1) = strcat(Ref{ii},' all');
            Ridx(end+1,1) = strcat(Ref{ii},' sst');
        end
    end   
    end
end



end
