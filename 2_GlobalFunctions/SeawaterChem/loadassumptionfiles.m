function assumptionfiles = loadassumptionfiles(Corrections)

% DEPRECATED pH -- using new approach of estimating local pH from ensemble
% if Corrections.phcorrection
%     load('PhanerozoicpHv5.mat')
%     assumptionfiles.GlobalpH = median(PhanerozoicpH,2);
% else
%     assumptionfiles.GlobalpH = zeros(101,1);
% end

% secular d18O 
load('Phanerozoicd18Ov5.mat')
assumptionfiles.Globald18O = median(PhanGlobalSW.GlobalSW,2);

% snowball earth
if Corrections.snowballearth
    assumptionfiles.Snowball = mean(PhanGlobalSW.Snowball,2);
else
    assumptionfiles.Snowball = zeros(101,1);
end

% predict locald18O
assumptionfiles.predictsw = load('SeawaterMdl.mat');

end