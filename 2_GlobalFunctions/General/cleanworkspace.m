function cleanworkspace(vars, keep, str)


if keep == "upper"
    keepvars = cellfun(@(v)v(1),isstrprop(vars,'upper'));
elseif keep == "lower"
    keepvars = cellfun(@(v)v(1),isstrprop(vars,'lower'));
elseif keep == "contains"
    keepvars = contains(vars,str);
elseif keep == "~contains"
    keepvars = ~contains(vars,str);
end

vars = vars(keepvars);
evalin("base", ['clearvars -except' sprintf(' %s',vars{:})])


end