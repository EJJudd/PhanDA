function [varout, dims] = combineruns(varin,idx,catdim)
% Set seed for reproduceability
rng(1,'twister')

varin = varin(:,idx);
varout = cell(size(varin,1),1);
dims = cellfun(@(x) size(x, catdim), varin);
dimmax = max(dims,[],'all');
for ii = 1:size(varin,1)
    varoutii = cell(1, sum(idx));
    for kk = 1:sum(idx)
       varoutii{kk} = varin{ii, kk};
       dimkk = size(varin{ii, kk}, catdim);
        if dimkk < dimmax
            varoutii{kk} = cat(catdim, varoutii{kk}, datasample(varin{ii, kk}, dimmax - dimkk, catdim));
        end
    end
    varout{ii,1} = cat(catdim, varoutii{:});
end

end





% function varout = combineruns(varin,idx,catdim)
% 
% varin = varin(:,idx);
% varout = cell(size(varin,1),1);
% for ii = 1:size(varin,1)
%     dim1 = size(varin{ii,1},catdim);
%     dim2 = size(varin{ii,2},catdim);
%     if dim1 > dim2
%         varout{ii,1} = cat(catdim,varin{ii,1},datasample(varin{ii,2},dim1,catdim));
%     elseif dim1 < dim2
%         varout{ii,1} = cat(catdim,datasample(varin{ii,1},dim2,catdim),varin{ii,2});
%     elseif dim1 == dim2
%         varout{ii,1} = cat(catdim,varin{ii,1},varin{ii,2});
%     end
% end
% 
% end