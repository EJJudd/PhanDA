function [TWave, TWstd] = timeweight(values,time)

% Step 2: Calculate time differences
timediff = diff(time);

% Step 3: Calculate weighted values
weightedvalues = values .* timediff;

% Step 4: Calculate time-weighted average
TWave = sum(weightedvalues) / sum(timediff);

% Step 5: Calculate squared differences
squareddiffs = (values - TWave).^2;

% Step 6: Calculate weighted squared differences
weightedsquareddiffs = squareddiffs .* timediff;

% Step 7: Calculate time-weighted variance
timeweightedvariance = sum(weightedsquareddiffs) / sum(timediff);

% Step 8: Calculate time-weighted standard deviation
TWstd = sqrt(timeweightedvariance);

end