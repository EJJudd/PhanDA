function d18cmacroEst = d18cmacroPSM_forward(sst,d18sw)

%Grossman 2012
a = 0.094;
b = -4.54;
c = 13.7-sst;

d18cmacroEst = ((-b - sqrt(b^2 - 4*a*c))./(2*a)) + d18sw;

% quadratic approximation of the O’Neil et al. (1969) equation by Hays and Grossman (1991)  
% d18cmacroEst = ((4.36-sqrt(11.4736+0.48*sst))./.24) + d18sw;

end
