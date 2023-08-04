function GMST = latweightgmst(temp,meandim)

if nargin == 1; meandim = [1,2]; end

% Adjust size of latweight to scale to size of temp
load("HadCM3Coordinates.mat","Lat")

if isa(temp,'double')
    ndim = ndims(temp);
    lwsize = [1,size(temp,2:ndim)];
    latweight = repmat(cosd(Lat), lwsize);
    latweight(isnan(temp)) = NaN;
    tempweight = temp .* latweight;
    GMST = squeeze( sum (tempweight,meandim,'omitnan') ./ ...
                sum (latweight,meandim,'omitnan') );
            
elseif isa(temp,'cell')
    GMST = cell(size(temp));
    for ii = 1:numel(temp)
        ndim = ndims(temp{ii});
        lwsize = [1,size(temp{ii},2:ndim)];
        latweight = repmat(cosd(Lat), lwsize);
        tempweight = temp{ii} .* latweight;
        latweight(isnan(tempweight)) = NaN;
        GMST{ii} = squeeze( sum (tempweight,meandim,'omitnan') ./ ...
                    sum (latweight,meandim,'omitnan') );
    end
    
end

end