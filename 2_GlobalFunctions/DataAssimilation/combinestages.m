function var = combinestages(var,varname,Preferences,endsize,catdim)

% Function to revise variable to account for combined stages

if varname == "GTS"
    for ii = 1:numel(Preferences.combstages(1,:))
       var.LowerBoundary(Preferences.combstages(1,ii)) = ...
            var.LowerBoundary(Preferences.combstages(2,ii));
        var.Average(Preferences.combstages(1,ii)) = mean([...
            var.UpperBoundary(Preferences.combstages(1,ii)), ...
            var.LowerBoundary(Preferences.combstages(1,ii))]);
        var.Stage(Preferences.combstages(1,ii)) = strcat(...
            var.Stage(Preferences.combstages(1,ii)),'/',...
            var.Stage(Preferences.combstages(2,ii)));
    end

elseif isnumeric(var) && numel(size(var)) == 2 && varname == "IceLine"
    s = size(var);
    for ii = 1:s(1)
        if any(ii == Preferences.combstages(1,:))
            var(ii,:) = mean([var(ii,:);var(ii+1,:)]);
        end
    end
    
elseif isnumeric(var) && numel(size(var)) == 2 && varname ~= "IceLine"
    s = size(var);
    var = [var,NaN(s)];
    for ii = 1:s(1)
        if any(ii == Preferences.combstages(1,:))
            var(ii,s(2)+1:end) = var(ii+1,1:s(2));
        else
            var(ii,s(2)+1:end) = var(ii,1:s(2));
        end
    end

elseif iscell(var)
    for ii = 1:numel(Preferences.combstages(1,:))
         var{Preferences.combstages(1,ii)} = cat(catdim, ...
            var{Preferences.combstages(1,ii)}, ...
            var{Preferences.combstages(2,ii)});
    end
end

var(Preferences.combstages(2,:),:) = [];
var = var(1:endsize,:);  


end