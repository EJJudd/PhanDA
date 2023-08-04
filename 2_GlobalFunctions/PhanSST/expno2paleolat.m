function Data = expno2paleolat(Data, PaleoCoordinates, ExpNo)

PaleoLat = NaN(height(Data),1);
PaleoLon = NaN(height(Data),1);
Data = addvars(Data,PaleoLat,PaleoLon,'After','ModLon');

uniqueagecoor = unique( table(Data.Stage, Data.ModLat, Data.ModLon, ...
    'VariableNames', {'Stage','ModLat','ModLon'}), 'rows');
for ii = 1:height(uniqueagecoor)
    % Index the data
    dataidx = find(Data.Stage == string(uniqueagecoor.Stage(ii)) & ...
                   Data.ModLat == uniqueagecoor.ModLat(ii) & ...
                   Data.ModLon == uniqueagecoor.ModLon(ii));
    % Index the coordinate lookup table
    tableidx = find(round(PaleoCoordinates.T0_0Ma(:,1),2) == round(uniqueagecoor.ModLat(ii),2) & ...
                    round(PaleoCoordinates.T0_0Ma(:,2),2) == round(uniqueagecoor.ModLon(ii),2));
    if isempty(tableidx)
    tableidx = find(round(PaleoCoordinates.T0_0Ma(:,1),1) == round(uniqueagecoor.ModLat(ii),1) & ...
                    round(PaleoCoordinates.T0_0Ma(:,2),1) == round(uniqueagecoor.ModLon(ii),1));
    end
    
    % Index the experiment numbers and make labels for the lookup table
    plateage = ExpNo.PlateAge(ExpNo.Stage == string(uniqueagecoor.Stage(ii)));
    platelab = sprintf("T%d_%dMa",floor(plateage),mod(plateage,1)*10);
    Data.PaleoLat(dataidx) = mean(PaleoCoordinates.(platelab)(tableidx,1));
    Data.PaleoLon(dataidx) = mean(PaleoCoordinates.(platelab)(tableidx,2));

end

if any(isnan(Data.PaleoLat))
    warning('You have NaN paleo-coordinate values')
end

end