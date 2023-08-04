function d18vpdb = vsmow2vpdb(d18vsmow)

% Based on the relationship described in Coplen et al., 2002 (p. 36):
    % Coplen, Tyler B. Compilation of minimum and maximum isotope ratios
    % of selected elements in naturally occurring terrestrial materials and 
    % reagents. Vol. 1. No. 4222. US Department of the Interior, USGS, 2002.
% And as recommended by the IUPAC in Brand et al., 2014 (Section 3.6):
    % Brand, Willi A., et al. "Assessment of international reference 
    % materials for isotope-ratio analysis (IUPAC Technical Report)." 
    % Pure and Applied Chemistry 86.3 (2014): 425-467.

d18vpdb = 0.97001 * d18vsmow - 29.99;

end