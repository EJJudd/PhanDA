function assumptionfiles = loadassumptionfiles(Corrections)

% DEPRECATED pH -- using new approach of estimating local pH from ensemble
% if Corrections.phcorrection
%     load('PhanerozoicpHv5.mat')
%     assumptionfiles.GlobalpH = median(PhanerozoicpH,2);
% else
%     assumptionfiles.GlobalpH = zeros(101,1);
% end

% secular d18O 
load('Phanerozoicd18Ov6.mat')
assumptionfiles.Globald18O = median(PhanGlobalSW.GlobalSW,2);

% snowball earth & Veizer sw corrections
if Corrections.swcorrection
    assumptionfiles.Snowball = mean(PhanGlobalSW.Snowball,2);
    assumptionfiles.Veizer = PhanGlobalSW.Veizer;
else
    assumptionfiles.Snowball = zeros(101,1);
    assumptionfiles.Veizer = zeros(101,1);
end

% predict locald18O
assumptionfiles.predictsw = load('SeawaterMdl.mat');

end