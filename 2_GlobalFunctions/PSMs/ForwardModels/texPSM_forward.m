function [tex86] = texPSM_forward(lat,lon,t,runname,type,stol,...
    alpha_samples,beta_samples,tau2_samples,Data_Input)
% lat = latitude of site
% lon = longitude of site
% t = temperatures
% varargin:
% 1) enter 'SST' or 'subT' to select calibration (default is SST)
% 2) enter "standard" or "analog" to specify calibration mode (default is standard)
% 3) search tolerance in t units (required for analog mode)

%% deal with optional arguments
% ng=nargin;
% if ng==6
%     runname=varargin{1};
%     type=varargin{2};
%     stol=varargin{3};
% elseif ng==5
%     error('To use analog mode, enter a search tolerance in TEX units');
% elseif ng==4
%     runname=varargin{1};
%     type="standard";
% elseif ng==3
%     runname='SST';
%     type="standard";
% else
%     error('You entered too many or too few arguments');
% end

%% error checking
if ~strcmp(type,"standard") && ~strcmp(type,"analog")
    error('please enter "analog" to specify analog mode')
end

%ensure vector
t=t(:);
lat=lat(:);
lon=lon(:);

%% Load data files needed in the analysis:

    
% if strcmp(type,"standard")
%     load(['tex/', 'Output_SpatAg_', runname, '/params_standard'],...
%         'alpha_samples_comp','beta_samples_comp','tau2_samples','Locs_Comp');
% else
%     load(['tex/', 'Output_SpatAg_', runname, '/params_analog'],...
%         'alpha_samples','beta_samples','tau2_samples');
%     load(['Data_Input_SpatAg_', runname],'Data_Input');
% end

% grid spacing is hard-coded here, will never change:
grid_half_space=10;
%% Figure out the alpha and beta series to draw from.       
if strcmp(type,"standard")
    Nloc = length(t);
    inder_g = NaN(Nloc,1);
    for i=1:Nloc
    inder_g(i)=find(abs(Locs_Comp(:,1)-lon(i))<=grid_half_space & abs(Locs_Comp(:,2)-lat(i))<=grid_half_space);
    end
    % Extract the alpha, beta series for that location:
    alpha_samples=alpha_samples_comp(inder_g, :);
    beta_samples=beta_samples_comp(inder_g, :);
else %analog option
    %find cells with modern T obs within the tolerance: 
    %number of big grids:
    N_bg = length(Data_Input.Locs(:,1));
    %NEW: calculate mean SSTs across spatial grid cells
    spatialMean = NaN(N_bg,1);
    for i=1:N_bg
        spatialMean(i)=mean(Data_Input.Target_Stack(Data_Input.Inds_Stack==i));
    end
    %identify mean values within the tolerance
    inder_g = spatialMean >= (mean(t)-stol) & spatialMean <= (mean(t)+stol);
    alpha_samples=alpha_samples(inder_g, :);
    beta_samples=beta_samples(inder_g, :);
    %tile tau2 to match
    tau2_samples=repmat(tau2_samples,1,size(alpha_samples,1));
    %reshape
    alpha_samples=reshape(alpha_samples,1,size(alpha_samples,1)*size(alpha_samples,2));
    beta_samples=reshape(beta_samples,1,size(beta_samples,1)*size(beta_samples,2));
end


%% Predict TEX86 values
tex86 = normrnd(t .* beta_samples + alpha_samples,repmat(sqrt(tau2_samples),length(t),1));
%downsample to 2000 iterations for ease
iters = length(alpha_samples);
randind = randsample(iters,1000);
tex86 = tex86(:,randind);
%any tex values outside 0 to 1 are forced to be in that range.
tex86(tex86>1)=1;
tex86(tex86<0)=0;


end


















%% GRAVEYARD: FORMER VERSION
%% deal with optional arguments
% ng=nargin;
% if ng==6
%     runname=varargin{1};
%     type=varargin{2};
%     stol=varargin{3};
% elseif ng==5
%     error('To use analog mode, enter a search tolerance in TEX units');
% elseif ng==4
%     runname=varargin{1};
% elseif ng==3
%     runname='SST';
%     type="standard";
% else
%     error('You entered too many or too few arguments');
% end
% 
% %% error checking
% if ~strcmp(type,"standard") && ~strcmp(type,"analog")
%     error('please enter "analog" to specify analog mode')
% end
% 
% %ensure vector
% t=t(:);
% 
% %% Load data files needed in the analysis:
% Ntk=20000;
%     load(['ModelOutput/', 'Output_SpatAg_', runname, '/tau2_samples'],'tau2_samples');
% if strcmp(type,"standard")
%     load(['ModelOutput/', 'Output_SpatAg_', runname, '/alpha_samples_comp'],'alpha_samples_comp');
%     load(['ModelOutput/', 'Output_SpatAg_', runname, '/beta_samples_comp'],'beta_samples_comp');
%     load(['ModelOutput/', 'Output_SpatAg_', runname, '/Locs_Comp'],'Locs_Comp');
% else
%     load(['ModelOutput/', 'Output_SpatAg_', runname, '/alpha_samples'],'alpha_samples')
%     load(['ModelOutput/', 'Output_SpatAg_', runname, '/beta_samples'],'beta_samples')
%     load(['ModelOutput/Data_Input_SpatAg_', runname],'Data_Input');
%     alpha_samples=[alpha_samples.field];
%     beta_samples=[beta_samples.field];
%     %trim to 20000
%     alpha_samples=alpha_samples(:,end-Ntk+1:end);
%     beta_samples=beta_samples(:,end-Ntk+1:end);
% end
% 
% % grid spacing is hard-coded here, will never change:
% grid_half_space=10;
% % trim tau^2 as it may not have had burnin removed:
% tau2_samples=tau2_samples(end-Ntk+1:end)';
% %% Figure out the alpha and beta series to draw from.       
% if strcmp(type,"standard")
%     inder_g=find(abs(Locs_Comp(:,1)-lon)<=grid_half_space & abs(Locs_Comp(:,2)-lat)<=grid_half_space);
%     % Extract the alpha, beta series for that location:
%     alpha_samples = mean(alpha_samples_comp(inder_g, :),1);
%     beta_samples = mean(beta_samples_comp(inder_g, :),1);
% else %analog option
%     %cycle through the alpha/beta grid cells, find those that feature mean
%     %modern T obs within the tolerance: 
%     
%     %number of big grids:
%     N_bg=length(Data_Input.Locs(:,1));
%     inder_g=[];
% 
%     for kk=1:N_bg
%         %find the SST obs corresponding to this index location:
%         vals=Data_Input.Target_Stack(Data_Input.Inds_Stack==kk);
%         %if the mean of vals in the big grids within tolerance, add it to inder_g
%         if mean(vals)>=(mean(t)-stol) && mean(vals)<=(mean(t)+stol)
%             inder_g=[inder_g; kk];
%         end
%     end
%     %error check in case no analogs were found
%     if isempty(inder_g)
%         error('No analogs were found. Make your search tolerance wider.')
%     else
%     end
%     alpha_samples=alpha_samples(inder_g, :);
%     beta_samples=beta_samples(inder_g, :);
%     %tile tau2 to match
%     tau2_samples=repmat(tau2_samples,1,size(alpha_samples,1));
%     %reshape
%     alpha_samples=reshape(alpha_samples,1,size(alpha_samples,1)*size(alpha_samples,2));
%     beta_samples=reshape(beta_samples,1,size(beta_samples,1)*size(beta_samples,2));
% end
% 
% 
% %% Predict TEX86 values
% tex86 = normrnd(t * beta_samples + alpha_samples,repmat(sqrt(tau2_samples),length(t),1));
% %downsample to 2000 iterations for ease
% iters = length(alpha_samples);
% randind = randsample(iters,1000);
% tex86 = tex86(:,randind);
% %any tex values outside 0 to 1 are forced to be in that range.
% tex86(tex86>1)=1;
% tex86(tex86<0)=0;