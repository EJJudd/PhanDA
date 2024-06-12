%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Hothouse Habitability  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sup. Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mdldir = '/Users/emilyjudd/Library/Mobile Documents/com~apple~CloudDocs/ModelOutputs/ExpFiles/';
load("HadCM3Coordinates.mat")
t = 45;
% Induan
iTAS = kel2cel(squeeze(ncread([mdldir,'scotese_4co2a_059_252Ma_1120pCO2_M1_tas.nc'],'tas')));
iMASK = ncread([mdldir,'scotese_4co2a_059_252Ma_1120pCO2_M1_tos.nc'],'tos_anomaly');
iMASK(~isnan(iMASK)) = 1; iMASK(isnan(iMASK)) = 0;
iWMMT = max(iTAS,[],3);
iAT = iWMMT; iAT(iAT<t) = NaN;
iAT(iAT<50) = 45;
iAT(iAT>= 50) = 50;
% PETM
pTAS = kel2cel(squeeze(ncread([mdldir,'scotese_4co2a_099_52Ma_1120pCO2_M1_tas.nc'],'tas')));
pMASK = ncread([mdldir,'scotese_4co2a_099_52Ma_1120pCO2_M1_tos.nc'],'tos_anomaly');
pMASK(~isnan(pMASK)) = 1; pMASK(isnan(pMASK)) = 0;
pWMMT = max(pTAS,[],3);
pAT = pWMMT; pAT(pAT<t) = NaN;
pAT(pAT<50) = 45;
pAT(pAT>= 50) = 50;
% Turonian
tTAS = kel2cel(squeeze(ncread([mdldir,'scotese_4co2a_091_92Ma_1120pCO2_M1_tas.nc'],'tas')));
tMASK = ncread([mdldir,'scotese_4co2a_091_92Ma_1120pCO2_M1_tos.nc'],'tos_anomaly');
tMASK(~isnan(tMASK)) = 1; tMASK(isnan(tMASK)) = 0;
tWMMT = max(tTAS,[],3);
tAT = tWMMT; tAT(tAT<t) = NaN;
tAT(tAT<50) = 45;
tAT(tAT>= 50) = 50;

%% Make figure
cm = hex2rgb({'#ffb703';'#9b2226'},1);

figure('Position',[326 450 1102 278],'Color','w')
tiledlayout(1,3,'Padding','none','TileSpacing','compact');
nexttile
ax = worldmap('World');mlabel off, plabel off, gridm off
contourm(Lat,Lon,iMASK,.5,'k','LineWidth',1.5)
pcolorm(Lat,Lon,iAT)
colormap(cm)
caxis([45 50])
textm(0,-165,'WMMT ≥ 45','FontWeight','bold',...
   'FontName','helvetica','fontsize',11,'color',cm(1,:))
textm(-10,-155,sprintf('(N = %.0f)',numel(find(iAT>=40))),...
   'FontName','helvetica','fontsize',10,'color',cm(1,:))
textm(-25,-165,'WMMT ≥ 50','FontWeight','bold',...
   'FontName','helvetica','fontsize',11,'color',cm(2,:))
textm(-35,-155,sprintf('(N = %.0f)',numel(find(iAT>=48))),...
   'FontName','helvetica','fontsize',10,'color',cm(2,:))
title(['Induan (~251 Ma)',newline,sprintf('GMST = %.f',...
    latweightgmst(mean(iTAS,3))),char(176),'C'],...
    'FontWeight','bold','FontName','helvetica','fontsize',11)

nexttile
ax = worldmap('World');mlabel off, plabel off, gridm off
contourm(Lat,Lon,tMASK,.5,'k','LineWidth',1.5)
pcolorm(Lat,Lon,tAT)
caxis([45 50])
textm(0,-165,'WMMT ≥ 45','FontWeight','bold',...
   'FontName','helvetica','fontsize',11,'color',cm(1,:))
textm(-10,-155,sprintf('(N = %.0f)',numel(find(tAT>=40))),...
   'FontName','helvetica','fontsize',10,'color',cm(1,:))
textm(-25,-165,'WMMT ≥ 50','FontWeight','bold',...
   'FontName','helvetica','fontsize',11,'color',cm(2,:))
textm(-35,-155,sprintf('(N = %.0f)',numel(find(tAT>=48))),...
   'FontName','helvetica','fontsize',10,'color',cm(2,:))
title(['Turonian (~92 Ma)',newline,sprintf('GMST = %.f',...
    latweightgmst(mean(tTAS,3))),char(176),'C'],...
    'FontWeight','bold','FontName','helvetica','fontsize',11)

nexttile
ax = worldmap('World');mlabel off, plabel off, gridm off
contourm(Lat,Lon,pMASK,.5,'k','LineWidth',1.5)
pcolorm(Lat,Lon,pAT)
caxis([45 50])
textm(0,-165,'WMMT ≥ 45','FontWeight','bold',...
   'FontName','helvetica','fontsize',11,'color',cm(1,:))
textm(-10,-155,sprintf('(N = %.0f)',numel(find(pAT>=40))),...
   'FontName','helvetica','fontsize',10,'color',cm(1,:))
textm(-25,-165,'WMMT ≥ 50','FontWeight','bold',...
   'FontName','helvetica','fontsize',11,'color',cm(2,:))
textm(-35,-155,sprintf('(N = %.0f)',numel(find(pAT>=48))),...
   'FontName','helvetica','fontsize',10,'color',cm(2,:))
title(['PETM (~56 Ma)',newline,sprintf('GMST = %.f',...
    latweightgmst(mean(pTAS,3))),char(176),'C'],...
    'FontWeight','bold','FontName','helvetica','fontsize',11)
