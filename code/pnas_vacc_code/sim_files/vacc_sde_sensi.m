function [ data_out ] = vacc_sde_sensi( y0, tmax, k, choose_risk, wl, wh, d, num_realisations, sigma )
% -----------------------------------------
%
% MODIFICATION OF VACC_SDE
%
% Now parameters are sampled from a triangular distribution
% Risk function may choose between linear or pyramid
%
% -----------------------------------------------------

% Original parameter values for Measles - time unit years
u=1/50;                 % Birth - Death rate
R0=16;                  % R0 value
g=365/13;               % Recovery rate
a=0;                    % Optional immigration rate

% Create triangular distributions to sample parameters from
u_dist = makedist('Triangular','a',0.5*u,'b',u,'c',1.5*u);
R0_dist = makedist('Triangular','a',0.5*R0,'b',R0,'c',1.5*R0);
g_dist = makedist('Triangular','a',0.5*g,'b',g,'c',1.5*g);
% d_dist = makedist('Triangular','a',0.99*d,'b',d,'c',1.01*d);
% d changes location of critical transition
k_dist = makedist('Triangular','a',0.5*k,'b',k,'c',1.5*k);

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
        
        elseif choose_risk == 3 
        % Linear increase function
        wfun = @(t) (wh-wl)*t/tmax + wl;
        
    end
        
        % Discretize to values for plotting purposes
        wvals=zeros(size(t));
    for i=1:length(t)
        wvals(i)=wfun(t(i));
    end
    

% White noise
G = @(t,y) sigma*[1;1e-3;1];

% Initialise array for output data
% This will have columns [t,w,y1,y2,y3,...,yn]  (n realisations)
data_out=zeros([length(t),3*num_realisations+2]);
    % Time column
    data_out(:,1)=t;
    % Risk function column
    data_out(:,2)=wvals;
    
    
    
% Begin realisations
for realisation_count = 1:num_realisations
        
    % Sample parameters from the distributions
    u = u_dist.random(1);
    R0 = R0_dist.random(1);
    g = g_dist.random(1);
    b=R0*(g+u);
    
    % d = d_dist.random(1);
    k = k_dist.random(1);
    
    
    % Evolution equations with sampled parameter values
    system = @(t,y)[u*(1-y(3))-b*y(1)*y(2)-u*y(1);...
                b*y(1)*y(2)-g*y(2)-u*y(2)+a;...
                k*y(3)*(1-y(3))*(-wfun(t)+y(2)+d*(2*y(3)-1))-a];

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

