function [GMST, LTG, LTGsst, TASpost, TOSpost, idx, Cgrid, Ycov] = runkf(...
    ensusemeta, ensusemat, Yuse, Ruse, Yeuse, idx, assvar, covanalysis)

dimorder = ["lat";"lon"];
TASpost = [];
TOSpost = [];
GMST = [];
Cgrid = cell(size(idx,1),1);
Ycov = cell(size(idx,1),1);
% (c) Assimilate each idx configuration
for ii = 1:size(idx,1)
    kf = kalmanFilter;
    kf = kf.prior(ensusemat(:,idx{ii,1}));
    kf = kf.observations(Yuse);
    kf = kf.uncertainties(Ruse);
    kf = kf.estimates(Yeuse(:,idx{ii,1}));
    kf = kf.deviations("return");
    output = kf.run;
    tasmean = ensusemeta.regrid("tas", output.Amean, 'order', dimorder);
    tasdev = ensusemeta.regrid("tas", output.Adev, 'order', dimorder);
    taspost = kel2cel(tasmean+tasdev);
    ltg = squeeze(mean(taspost,2));
    idx{ii,3} = thresholdltg(ltg);
    TASpost = cat(3,TASpost,taspost);
    GMST = cat(1, GMST, latweightgmst(taspost));
    
    if any(contains(assvar,"tos"))
        tosmean = ensusemeta.regrid("tos", output.Amean, 'order', dimorder);
        tosdev = ensusemeta.regrid("tos", output.Adev, 'order', dimorder);
        tospost = tosmean+tosdev;
        tospost(tospost<-1.8) = -1.8;
        TOSpost = cat(3,TOSpost,tospost);
    end
    
    if covanalysis
       [C, Ycov{ii}] = kf.covariance;
        Cgrid{ii} = ensusemeta.regrid("tas", C, 'order', dimorder);
    end
    
end

LTG = squeeze(mean(TASpost,2));
LTGsst = squeeze(nanmean(TOSpost,2));

end




