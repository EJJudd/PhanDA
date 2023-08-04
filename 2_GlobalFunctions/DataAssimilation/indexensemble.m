function [idx, exprun, exps, labs] = indexensemble(ensusemeta,a)

% Load experiment file
load("ExpNo.mat")

% Define experiment names and numbers
dims = string(ensusemeta.ensembleDimensions{1});
exprunall = ensusemeta.ensemble{1,1}{1, dims == "exprun"};
runs = unique(exprunall);
newruns = ["scotese06","scotese07","scotese08"];
orgruns = runs(~contains(runs,newruns));
expno = ensusemeta.ensemble{1,1}{1, dims == "expno"};
no1 = ExpNo.ExpNo1(a);

% index those experiments
exprun = cell(35,1);
idx = cell(35,2);
exps = cell(35,1);
labs = cell(35,1);

c = 1;
% 1st configuration: all runs, both time steps
exprun{c} =  exprunall(contains(exprunall,runs));
idx{c,1} = find(contains(exprunall,runs));
idx{c,2} = [1:numel(idx{c,1})]';
labs{c} = sprintf("Time steps: 2%sPriors: All",newline);
c=c+1;
% 2nd configuration: original runs, both time steps
exprun{c} = exprunall(contains(exprunall,orgruns));
idx{c,1} = find(contains(exprunall,orgruns));
idx{c,2} = [max(idx{c-1,2})+1:max(idx{c-1,2})+numel(idx{c,1})]';
labs{c} = sprintf("Time steps: 2%sPriors: Orig.",newline);
c=c+1;
% 3rd configuration: new runs, both time steps
exprun{c} = exprunall(contains(exprunall,newruns));
idx{c,1} = find(contains(exprunall,newruns));
idx{c,2} = [max(idx{c-1,2})+1:max(idx{c-1,2})+numel(idx{c,1})]';
labs{c} = sprintf("Time steps: 2%sPriors: New",newline);
c=c+1;
% 4th configuration: all runs, best time step
exprun{c} = exprunall(contains(exprunall,runs)&expno==no1);
idx{c,1} = find(contains(exprunall,runs)&expno==no1);
idx{c,2} = [max(idx{c-1,2})+1:max(idx{c-1,2})+numel(idx{c,1})]';
labs{c} = sprintf("Time steps: 1%sPriors: All",newline);
c=c+1;
% 5th configuration: original runs, best time step
exprun{c} = exprunall(contains(exprunall,orgruns)&expno==no1);
idx{c,1} = find(contains(exprunall,orgruns)&expno==no1);
idx{c,2} = [max(idx{c-1,2})+1:max(idx{c-1,2})+numel(idx{c,1})]';
labs{c} = sprintf("Time steps: 1%sPriors: Orig",newline);
c=c+1;
% 6th configuration: new runs, best time step
exprun{c} = exprunall(contains(exprunall,newruns)&expno==no1);
idx{c,1} = find(contains(exprunall,newruns)&expno==no1);
idx{c,2} = [max(idx{c-1,2})+1:max(idx{c-1,2})+numel(idx{c,1})]';
labs{c} = sprintf("Time steps: 1%sPriors: New",newline);
c=c+1;
% 8: cycling through all runs, leaving one out (both time steps)
for ii = 1:numel(runs)
    str = runs; lo = str(ii); str(ii) = [];
    exprun{c} = exprunall(contains(exprunall,str));
    idx{c,1} = find(contains(exprunall,str));
    idx{c,2} = [max(idx{c-1,2})+1:max(idx{c-1,2})+numel(idx{c,1})]';
    labs{c} = sprintf("Time steps: 2%sPriors: All (no %s)",newline,lo);
    c=c+1;
end
% 8: cycling through all runs, leaving one out (best time step)
for ii = 1:numel(runs)
    str = runs; lo = str(ii); str(ii) = [];
    exprun{c} = exprunall(contains(exprunall,str)&expno==no1);
    idx{c} = find(contains(exprunall,str)&expno==no1);
    idx{c,2} = [max(idx{c-1,2})+1:max(idx{c-1,2})+numel(idx{c,1})]';
    labs{c} = sprintf("Time steps: 1%sPriors: All (no %s)",newline,lo);
    c=c+1; 
end
% 5: cycling through original runs, leaving one out (both time steps)
for ii = 1:numel(orgruns)
    str = orgruns; lo = str(ii); str(ii) = [];
    exprun{c} = exprunall(contains(exprunall,str));
    idx{c} = find(contains(exprunall,str));
    idx{c,2} = [max(idx{c-1,2})+1:max(idx{c-1,2})+numel(idx{c,1})]';
    labs{c} = sprintf("Time steps: 2%sPriors: Orig. (no %s)",newline,lo);
    c=c+1;
end
% 5: cycling through original runs, leaving one out (best time step)
for ii = 1:numel(orgruns)
    str = orgruns; lo = str(ii); str(ii) = [];
    exprun{c} = exprunall(contains(exprunall,str)&expno==no1);
    idx{c} = find(contains(exprunall,str)&expno==no1);
    idx{c,2} = [max(idx{c-1,2})+1:max(idx{c-1,2})+numel(idx{c,1})]';
    labs{c} = sprintf("Time steps: 1%sPriors: Orig. (no %s)",newline,lo);
    c=c+1;    
end

% 3: cycling through new runs, leaving one out (both time steps)
for ii = 1:numel(newruns)
    str = newruns; lo = str(ii); str(ii) = [];
    exprun{c} = exprunall(contains(exprunall,str));
    idx{c} = find(contains(exprunall,str));
    idx{c,2} = [max(idx{c-1,2})+1:max(idx{c-1,2})+numel(idx{c,1})]';
    labs{c} = sprintf("Time steps: 2%sPriors: New (no %s)",newline,lo);
    c=c+1;
end

for ii = 1:numel(idx)/2
    exps{ii} = unique(exprun{ii});
end
        


end