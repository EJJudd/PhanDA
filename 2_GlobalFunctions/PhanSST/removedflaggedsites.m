function Data = removedflaggedsites(Data)

% This script removes sites that have been deemed unreliable/suspect based
% on manual inspection. A full explanation of why data from each site are
% ommitted can be found in the file "PhanDA_FlaggedSites.xlsx" 
% (PhanDA > 4_NonGlobalFiles > PhanSST)


% --> (1) Azmy and Jin, 2018
%         Sheguindah Shale and Stony Mountain Formation
Data.DiagenesisFlag(Data.Formation=="Stony Mountain" & Data.LeadAuthor=="Azmy") = 1;
Data.DiagenesisFlag(Data.Formation=="Sheguindah Shale" & Data.LeadAuthor=="Azmy") = 1;

% --> (2) Bruckschen et al., 2001 & Brand and Bruckschen, 2002
%         Serpukhovian Askyn data 
Data.Environment(Data.SiteName=="Askyn"&Data.Stage=="Serpukhovian") = "restricted";

% --> (3) Buggisch et al., 2008
%         Combine data from Tellego & Mirador de Vegamian
Data.PaleoLat(Data.SiteName=="Tellego"&Data.Stage=="Serpukhovian") = ...
    unique(Data.PaleoLat(Data.SiteName=="Mirador De Vegamian"&Data.Stage=="Serpukhovian"));
Data.PaleoLon(Data.SiteName=="Tellego"&Data.Stage=="Serpukhovian") = ...
    unique(Data.PaleoLon(Data.SiteName=="Mirador De Vegamian"&Data.Stage=="Serpukhovian"));

% --> (4) Frank et al., 2015
%         deglacial/interglacial data
Data.DiagenesisFlag(Data.PublicationDOI=="https://doi.org/10.1016/j.palaeo.2014.11.016" & ...
    (Data.Stage=="Sakmarian" | Data.Stage=="Artinskian") & Data.ProxyValue<-1.5) = 1;

% --> (5) Garbelli et al., 2016
%         Gyanyima data
Data.PaleoLat(Data.SiteName == "Gyanyima") = -40;
Data.PaleoLon(Data.SiteName == "Gyanyima") = 63.75;

% --> (6) Grossman et al., 2008 
%         Combine data from Fayetteville Fm.
Data.PaleoLat(Data.Formation=="Fayetteville") = -12.5;
Data.PaleoLon(Data.Formation=="Fayetteville") = -33.75;

% --> (7) Hornung et al., 2007 
%         Lanadian, early Carnian, and early Norian
Data.Environment(contains(Data.Stage,"Ladinian")&Data.LeadAuthor=="Hornung") = "restricted";
Data.Environment(Data.Stage=="Carnian"&Data.StagePosition=="early"&Data.LeadAuthor=="Hornung") = "restricted";
Data.Environment(Data.Stage=="Norian"&Data.StagePosition=="early"&Data.LeadAuthor=="Hornung") = "restricted";


% --> (8) Khelif et al., 2009 & 2014
%         Zanclean ODP978 data
Data.DiagenesisFlag(Data.SiteName=="ODP978"&Data.Stage=="Zanclean") = 1;


% --> (9) Le Houedec et al., 2013 
%         M'rirt data
Data.DiagenesisFlag(Data.SiteName=="M'rirt") = 1;

% --> (10) Mii et al., 1997 
%          data from the Kapp Starostin formation
Data.DiagenesisFlag(Data.Formation=="Kapp Starostin") = 1;

% --> (11) Tripati & Zachos, 2002
%          Data from Panama
Data.DiagenesisFlag(Data.LeadAuthor == "Tripati" & Data.Year==2002 & ...
    Data.Country == "Panama") = 1;

% --> (12) Wang et al., 2020
%          Phosphate data from the Xikou formation
Data.DiagenesisFlag(Data.SiteName=="Xikou" & Data.ProxyType=="d18p") = 1;

end
