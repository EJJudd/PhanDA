function [coasts, lat, lon, dem] = scotesecontourlines(mapfilenames,GTS)
    
figure, hold on

for jj = 1:numel(mapfilenames)
    stageName = string(GTS.Stage(jj));
    stageLab = sprintf('S%02d_%s',jj,stageName); 
    stageLab = strrep(stageLab,' ','');
    mapfilename = string(mapfilenames(jj));
    dem.(stageLab) = csvread(mapfilename,1,0);
    if length(dem.(stageLab)) == 65341
        lat.(stageLab) = reshape(dem.(stageLab)(:,2),361,181)';
        lon.(stageLab) = reshape(dem.(stageLab)(:,1),361,181)';
        dem.(stageLab) = reshape(dem.(stageLab)(:,3),361,181)';
    elseif length(dem.(stageLab)) == 64800
        lat.(stageLab) = reshape(dem.(stageLab)(:,2),360,180)';
        lon.(stageLab) = reshape(dem.(stageLab)(:,1),360,180)';
        dem.(stageLab) = reshape(dem.(stageLab)(:,3),360,180)';
    else
        disp('Error in dem size')
    end
    coasts.(stageLab) = contourm(lat.(stageLab),lon.(stageLab),dem.(stageLab),0,'k')';
    clf
end

close
    
end