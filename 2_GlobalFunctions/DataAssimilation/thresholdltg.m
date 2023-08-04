function   threshold = thresholdltg(ltg)

    load("HadCM3Coordinates.mat","Lat")
    x = mean(ltg,2);
    
    dt = .25 * range(x);
    
    % Criteria #1a
    % Temperatures between 80-90 N are warmer than temps between 60-70 N
    % plus .25 the total range of the LTG
    if mean(x(Lat>80&Lat<90)) > mean(x(Lat>60&Lat<70)) + dt 
        threshold = false;
    % Criteria #1b
    % Temperatures between 80-90 S are warmer than temps between 60-70 S
    % plus .25 the total range of the LTG
    elseif mean(x(Lat<-80&Lat>-90)) > mean(x(Lat<-60&Lat>-70)) + dt
        threshold = false;
    % Criteria #2
    % The change in temperature between sequential latitudes is >25% the
    % total range of the gradient 
    elseif any(abs(diff(x))/range(x)*100 > 33)
        threshold = false;
    % Criteria #3
    % Any polar temperatures > any tropical temperatures
    elseif max(x(abs(Lat)>66.5)) > min(x(abs(Lat)<23.5))
        threshold = false;
    else
        threshold = true;
    end
    
end
