function [ data_out ] = vacc_sde( y0, tmax, k, choose_risk, wl, wh, d, num_realisations, sigma )
% vacc_sde runs a stochastic simulation of the vaccine imitation model
%   data_out = vacc_sde( y0, tmax, k, choose_risk, wh, wl, d, num_realisations, noise )
%   Input arguments: 
%       y0 - intitial conition for y = [S,I,x]
%       infections
%       tmax - final time point of simulation
%       k - product of imitaiton rate and sensitivity to disease prevalence
%       choose_risk - type of risk function ( in {1,2,3,4,5} )
%       wh - maximum value of risk function
%       wl - minimum value of risk function
%       d - social pressure term
%       num_realisations - # realisations
%       sigma - strength of noise component ( in [0,1] )
%
%   Output variable:
%       data_out - is a time series represented as a matrix with cols
%       [t,w,y1,y2,...,yn] where yi is the ith realization

% -----------------------------------------------------

% Fixed parameter values for Measles - time unit years
u=1/50;                 % Birth - Death rate
R0=16;                  % R0 value
g=365/13;               % Recovery rate
b=R0*(g+u);             % Transmission rate
a=0;                    % Optional immigration rate

% Create time vector           
dt = 5e-4;
t = 0:dt:tmax;


% Assign risk function
    if choose_risk == 1
        % Pyramid function
        tspread=tmax/10;
        tcenter=tmax/2;
        wfun = @(t) (t<tcenter-tspread)*wl + ...
            (t>=tcenter-tspread && t<tcenter)*((wh-wl)*t/tspread+wh-tcenter*(wh-wl)/tspread)+...
            (t>=tcenter && t<tcenter+tspread)*(-(wh-wl)*t/tspread+wh+tcenter*(wh-wl)/tspread) +...
            (t>=tcenter+tspread)*wl;
       
    elseif choose_risk == 2
        % Box function
        ton = 100;
        toff = 200;
        wfun = @(t) (t<ton)*wl + (t>=ton && t<toff)*wh + (t>=toff)*wl;
           
    elseif choose_risk == 3 
        % Linear increase function
        wfun = @(t) (wh-wl)*t/tmax + wl;
        
    elseif choose_risk == 4
        % Flat then linear increase then flat
        tinc=tmax/8;
        wfun=@(t) (t<(tmax-tinc)/2).*wl + ...
            (t>=(tmax-tinc)/2 && t<(tmax+tinc)/2)*(wl+(t-(tmax-tinc)/2)*(wh-wl)/tinc) +...
            (t>=(tmax+tinc)/2)*wh;
        
    elseif choose_risk == 5
        % Gaussian curve
        c=100/tmax^2;
        wfun = @(t) (wh-wl)*exp(-c*(t-tmax/2)^2) + wl;
        
    elseif choose_risk == 6
        % Trapezoid
        tinc = tmax/10;
        ton = tmax/5;
        t1 = (tmax-ton)/2-tinc;
        t2 = t1+tinc;
        t3 = t2+ton;
        t4 = t3+tinc;
        wfun=@(t) (t<t1)*wl + ...
            (t>=t1 && t<t2)*(wl+(t-t1)*(wh-wl)/tinc) +...
            (t>=t2 && t<t3)* wh +...
            (t>=t3 && t<t4)*(wh-(t-t3)*(wh-wl)/tinc) +...
            (t>=t4)* wl;
    end
    
% Discretize to values for plotting purposes
    wvals=zeros(size(t));
for i=1:length(t)
    wvals(i)=wfun(t(i));
end


% Define the evolution equations (system of odes)
system = @(t,y)[u*(1-y(3))-b*y(1)*y(2)-u*y(1);...
                b*y(1)*y(2)-g*y(2)-u*y(2)+a;...
                k*y(3)*(1-y(3))*(-wfun(t)+y(2)+d*(2*y(3)-1))-a];


% White noise
G = @(t,y) sigma*[1;1e-3;1];

    
% Initialise array for output data
% This will have columns [t,w,y1,y2,y3,...,yn]  (n realisations)
% Each yi has form [S,I,x]

data_out=zeros([length(t),3*num_realisations+2]);
    % Time column
    data_out(:,1)=t;
    % Risk function column
    data_out(:,2)=wvals;
    
    
    
% Begin realisations
for realisation_count = 1:num_realisations

    % Random number generator seed
    seed = realisation_count; 
    
    % Implement sdetools
    ops = sdeset('SDEType','Ito','RandSeed',seed, 'ConstGFUN','yes','NonNegative','yes','UBoundOne','yes');
    y = sde_euler(system,G,t,y0,ops);

    % Save current data to the output array
    data_out(:,3*realisation_count : 3*realisation_count+2)=y;
    
    % Print complete
    confirm=['Realisation ',num2str(realisation_count),' complete'];
    disp(confirm)
end      


end

