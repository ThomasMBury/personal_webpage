function [ data_filt ] = time_filter( data, tinit, tfin, num_comps )
% time_filter extracts equally spaced points from a given time series
%
%   out = time_filter( data, tinit, tfin, num_comps )
%
%   Input arguments:
%       data - time-series in form of matrix with cols [t,y1,y2,...,yn]
%       tinit - desired intial time
%       tfin - desired final time
%       num_comps - desired # components for out
%
%   Output variable:
%       data_filt - filtered time series with cols [t,w,y1,y2,...,yn]



% Number of time components in original data
orig_num_comps = size(data,1);

% Largest time in original data
tmax = data(end,1);

% Original spacing
dt = tmax/(orig_num_comps+1);

% Desired range
range = tfin-tinit;


% Required spacing (number of comps to skip) between each recorded time in order to get desired
% number of components
spacing = floor((orig_num_comps/num_comps)*(range/tmax));
    
% Create filtered array
data_filt=zeros([num_comps,size(data,2)]);
for i=1:num_comps
        data_filt(i,:)=data(floor(tinit/dt+1) + (i-1)*spacing,:);
end

end

