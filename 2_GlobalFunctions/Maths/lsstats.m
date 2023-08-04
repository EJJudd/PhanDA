function LS = lsstats(lsmask)

load("HadCM3Coordinates.mat")
LatWeight = cosd(Lat);
s = size(lsmask,1);
% Preallocate structure
LS.lssum = NaN(numel(Lat),s);
LS.lssum_scaled = NaN(numel(Lat),s);
LS.lssumweight = NaN(numel(Lat),s);
LS.lssumweight_scaled = NaN(numel(Lat),s);

for a = 1:s
    lsmask{a} = mean(lsmask{a},3);
    lsmask{a}(1,:) = lsmask{a}(2,:);
    lsmask{a}(end,:) = lsmask{a}(end-1,:);
    LS.lssum(:,a) = sum(lsmask{a},2);
    LS.lssum_scaled(:,a) = 100*sum(lsmask{a},2)./sum(lsmask{a}(:));
    LS.lssumweight(:,a) = sum(lsmask{a},2).*LatWeight;
    LS.lssumweight_scaled(:,a) = 100*sum(lsmask{a},2).*LatWeight./(sum(sum(lsmask{a},2).*LatWeight));
end

end