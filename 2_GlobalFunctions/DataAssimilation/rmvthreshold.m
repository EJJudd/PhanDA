function [gmstfinal,ltgfinal,ltgsstfinal,taspostfinal] = ...
    rmvthreshold(gmst,ltg,ltgsst,taspost,idx)

passed = vertcat(idx{cell2mat(idx(:,3)) == 1,2});
gmstfinal = gmst(passed,1);
ltgfinal = ltg(:,passed);
ltgsstfinal = ltgsst(:,passed);
taspostfinal = taspost(:,:,passed);

end