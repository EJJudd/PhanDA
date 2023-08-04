function [omega,ph] = omgph(lat,lon,depth)
% function [omega,ph] = omgph(lat,lon,depth)
%
% This function gets calcite omega for a given core site(s) (at core
% site depth) as well as surface (0 m water depth) pH for use with the
% BAYMAG forward and prediction models. It draws these values from GLODAPv2
% (2016) unless the site is in the Gulf of Mexico, in which case the values
% are derived from GOMECC-2 observations.
%
% ----- Input -----
% lat = latitude in decimal degrees (scalar or N x 1 vector)
% lon = longitude in decimal degrees (scalar or N x 1 vector, -180 to 180)
% depth = water depth of core site in meters as a POSITIVE number (scalar
% or vector)
%
% ----- Output -----
% omega = omega (saturation) for calcite (scalar or vector)
% ph = surface water pH value (scalar or vector)
%
% ----- Dependencies -----
% 1) the GLODAPv2 netcdfs for OmegaC and pHtsinsitutp
% 2) EarthChordDistances_2.m
% 3) gom.mat - .mat file containing GOMECC-2 values
% 
% function created by Dr. Jessica Tierney, The University of Arizona (2019)

%% set-up
%get Nobs
Nobs=length(lat);
%put together entered locations
coords=[lon lat];

%load GLODAPv2 Omega and Surface pH
lat_f=ncread('GLODAPv2.2016b.OmegaC.nc','lat');
lon_f=ncread('GLODAPv2.2016b.OmegaC.nc','lon');
wdepth=ncread('GLODAPv2.2016b.OmegaC.nc','Depth');
%load whole field
omega_field=ncread('GLODAPv2.2016b.OmegaC.nc','OmegaC');
%only load the top layer
ph_field=ncread('GLODAPv2.2016b.pHtsinsitutp.nc','pHtsinsitutp',[1 1 1],[Inf Inf 1]);
%reorganize longitude and fields
Nlon=length(lon_f);
lon_f=[lon_f(Nlon/2-20+1:1:(Nlon-20))-360; lon_f((Nlon-20)+1:1:Nlon)-360; lon_f(1:1:Nlon/2-20)];
ph_field=[ph_field(Nlon/2-20+1:1:Nlon-20,:); ph_field(Nlon-20+1:1:Nlon,:); ph_field(1:1:Nlon/2-20,:)];
omega_field=[omega_field(Nlon/2-20+1:1:Nlon-20,:,:); omega_field(Nlon-20+1:1:Nlon,:,:); omega_field(1:1:Nlon/2-20,:,:)];

%manually cut out westernmost caribbean data to ensure sites to the west of
%it get assigned appropriately:
ph_carib=ph_field(108,103:108);
lat_carib=lat_f(103:108);
omega_carib=squeeze(omega_field(108,103:108,:));

%vectorize locations
[A,B] = meshgrid(lon_f,lat_f);
c=cat(2,A',B');
locs=reshape(c,[],2);

%reorganize into vector format
n_lon=length(lon_f);
n_lat=length(lat_f);
n_depth=length(wdepth);

omega_vec=reshape(omega_field,n_lon*n_lat,n_depth);
ph_vec=reshape(ph_field,n_lon*n_lat,1);

%pull out locations, data where there are obs for pH
locs_obs_ph=locs(~isnan(ph_vec),:);
ph_obs=ph_vec(~isnan(ph_vec));

%pull out locations, data where there is omega and put into a 33 x 1
%structure:
locs_obs_omega=cell(33,1);
omega_obs=cell(33,1);
for i=1:length(wdepth)
locs_obs_omega{i}=locs(~isnan(omega_vec(:,i)),:);
omega_obs{i}=omega_vec(~isnan(omega_vec(:,i)),i);
end

%load Gulf of Mexico data
gom=load('gom.mat');
%set poly for Gulf of Mexico, flag locations in it
poly_gom=[...
    -96.5 16.5
    -100.3 30.5
    -82 30.5
    -80.5 23];
ind_g=inpolygon(lon,lat,poly_gom(:,1),poly_gom(:,2));
%set poly for Caribbean Sea, flag locations in it
poly_carib=[...
    -77.5 8
    -90.8 18.6
    -82.4 22.9
    -72.5 17.5
    -72.5 8.8];
ind_c=inpolygon(lon,lat,poly_carib(:,1),poly_carib(:,2));
%% cycle through each lat, lon, depth and extract values
%preallocate vectors for omega and pH
omega=NaN(Nobs,1);
ph=NaN(Nobs,1);
dists_ph=NaN(Nobs,1);
dists_om=NaN(Nobs,1);
d_diff=NaN(Nobs,1);
%use a cutoff distance, esp. important for omega.
max_dist=700;

%now loop across coords to get data.
for i=1:Nobs
    tloc=coords(i,:);
    %first check for GoM locations
    if ind_g(i)
        d1=findnearest(depth(i),gom.depth);
        omega(i)=gom.omega(d1(1));
        %notice we are grabbing surface ph
        ph(i)=gom.ph(1);
        dists_om(i)=0;
        dists_ph(i)=0;
    %next check for Carib location
    elseif ind_c(i)
        d1=findnearest(depth(i),wdepth);
        d1=d1(1);
        lat1=findnearest(tloc(2),lat_carib);
        ph(i)=ph_carib(lat1(1));
        omega(i)=omega_carib(lat1(1),d1);
        %loop in case you grab a NaN for omega
        while isnan(omega(i))
            d1=d1-1;
            omega(i)=omega_carib(lat1(1),d1);
        end
        dists_om(i)=0;
        dists_ph(i)=0;    
    else
        %find ph, no cutoff
        [dmin,imin]=min(EarthChordDistances_2(locs_obs_ph,tloc));
        ph(i)=ph_obs(imin);
        %record distance from location
        dists_ph(i)=dmin;
        %find omega - first check for NaN depth
        if isnan(depth(i))
            omega(i)=NaN;
        else
        dnew=findnearest(depth(i),wdepth);
        dnew=dnew(1);
        [dmin,imin]=min(EarthChordDistances_2(locs_obs_omega{dnew},tloc));
        %if omega is greater than cutoff, search nearby depths
        while dmin > max_dist
            dnew=dnew-1;
            [dmin,imin]=min(EarthChordDistances_2(locs_obs_omega{dnew},tloc));
        end
        %record distances from location (depth and lat/lon), if you want 
        %you can output these
        d_diff(i)=depth(i)-wdepth(dnew);
        dists_om(i)=dmin;
        omega_now=omega_obs{dnew};
        omega(i)=omega_now(imin);
        end
    end
end