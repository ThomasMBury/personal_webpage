function [ data_out ] = vacc_sde_erlang( y0, tmax, k, choose_risk, wl, wh, d, num_realisations, sigma )
% vacc_sde runs a stochastic simulation of the vaccine imitation model
% WITH Erlang distributed infectious period
%
%   data_out = vacc_sde( y0, tmax, k, choose_risk, wh, wl, d, num_realisations, sigma )
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
n=20;                   % Erlang distribution shape parameter (number of I compartments)


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


% Define the evolution equations (system of odes) - just went crazy trying
% to do this for general n - so now putting in directly - i'm going to hell
% for this
system = @(t,y)[u*(1-y(22))-b*y(1)*sum(y(2:21))-u*y(1);...  %y(1)
                b*y(1)*sum(y(2:21))-(n*g + u)*y(2);...      %y(2)
                n*g*y(2)- (n*g + u)*y(3);...        %y(3)
                n*g*y(3)- (n*g + u)*y(4);...        %y(4)
                n*g*y(4)- (n*g + u)*y(5);...        %y(5)
                n*g*y(5)- (n*g + u)*y(6);...        %y(6)
                n*g*y(6)- (n*g + u)*y(7);...        %y(7)
                n*g*y(7)- (n*g + u)*y(8);...        %y(8)
                n*g*y(8)- (n*g + u)*y(9);...        %y(9)
                n*g*y(9)- (n*g + u)*y(10);...       %y(10)
                n*g*y(10)- (n*g + u)*y(11);...      %y(11)
                n*g*y(11)- (n*g + u)*y(12);...      %y(12)
                n*g*y(12)- (n*g + u)*y(13);...      %y(13)
                n*g*y(13)- (n*g + u)*y(14);...      %y(14)
                n*g*y(14)- (n*g + u)*y(15);...      %y(15)
                n*g*y(15)- (n*g + u)*y(16);...      %y(16)
                n*g*y(16)- (n*g + u)*y(17);...      %y(17)
                n*g*y(17)- (n*g + u)*y(18);...      %y(18)
                n*g*y(18)- (n*g + u)*y(19);...      %y(19)
                n*g*y(19)- (n*g + u)*y(20);...      %y(20)
                n*g*y(20)- (n*g + u)*y(21);...      %y(21)
                k*y(22)*(1-y(22))*(-wfun(t)+sum(y(2:21))+d*(2*y(22)-1))];%y(22)


% White noise
noise_vec = zeros(n+2,1);
noise_vec(1) = 1; % S noise
noise_vec(n+2) = 1; % x noise
noise_vec(2:n+1) = ones(n,1)*1e-3;
G = @(t,y) sigma*noise_vec;

    
% Initialise array for output data
% This will have columns [t,w,y1,y2,y3,...,yn]  (n realisations)
% Each yi has form [S,I,x]

data_out=zeros([length(t),22*num_realisations+2]);
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
    data_out(:,22*realisation_count-19 :22*realisation_count+2)=y;
    
    % Print complete
    confirm=['Realisation ',num2str(realisation_count),' complete'];
    disp(confirm)
end      


end



