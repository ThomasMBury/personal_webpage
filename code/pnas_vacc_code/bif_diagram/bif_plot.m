% Script to make bifurcation plot for vaccine model
close all
clear


%---------------------
% Set up plot
%---------------------

% Properties
dw = 1e-6;      % Resolution
w_fin = 7e-4;   % Omega range
lw = 1.5;         % Line width of plot


% Fixed parameter values for Measles - time unit years
u=1/50;                 % Birth - Death rate
R0=16;                  % R0 value
g=365/13;               % Recovery rate
b=R0*(g+u);             % Transmission rate

% Social parameters
d = 5e-4;

% Intersection points on bif plot
w_int1 = d;
w_int2 = d*(1-2/R0);
w_int3 = (u/b)*(R0-1) - d;


% Construct e1
w1f = 0:dw:w_int1;
w1s = w_int1:dw:w_fin;
e1f = w1f.^0;
e1s = w1s.^0;


% Construct e3
w3 = 0:dw:w_int1;
e3 = (1/2)*(1+w3/d);


% Construct e4
w4f = 0:dw:w_int3;
w4s = w_int3:dw:w_fin;
e4f = 0*w4f.^0;
e4s = 0*w4s.^0;


% Construct e5
w5 = w_int3:dw:w_int2;
e5 = (u*(1-1/R0) - (d+w5)*(u+g))/(u-2*d*(u+g));



%------------------------
% Make plot
%------------------------

figure()
hold on

% Set figure size
set(gcf,'Units','Inches')
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 10.5 2])

% Axes properties
axis([0,w_fin,-0.01,1])
set(gca,'FontSize',12)


ax = gca;
ax.XTick = [0,1e-4,2e-4,3e-4,4e-4,5e-4,6e-4];
ax.XTickLabel = {'0','1','2','3','4','5','6'};

xlabh = get(gca,'XLabel');
set(xlabh,'Position',[6.5e-4,-.05,0])

xlabel('$w \times 10^{-4}$','Interpreter','latex')
ylabel('Uptake')


% plot e1
plot(w1f,e1f,'Color','k','DisplayName','e1','LineWidth',lw)
plot(w1s,e1s,'--','Color','k','LineWidth',lw)


% plot e3
plot(w3,e3,'--','Color','k','DisplayName','e3','LineWidth',lw)


% plot e4
plot(w4f,e4f,'--','Color','k','DisplayName','e4','LineWidth',lw)
plot(w4s,e4s,'-','Color','k','LineWidth',lw)


% plot e5
plot(w5,e5,'--','Color','k','LineWidth',lw)










