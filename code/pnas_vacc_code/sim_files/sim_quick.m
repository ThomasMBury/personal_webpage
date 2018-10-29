% Run a single simulation

clear
close all

% Input variables - see vacc_sde.m for definitions
y0 = [0.01,1e-4,0.99];
tmax = 100;
k = 500;
choose_risk = 3;
wl = 1e-4;              
wh = 7e-4;
d = 5e-4;
num_realisations = 5;     
sigma = 0.01;     

% Run vacc_sde
% Output is matrix with cols [t,w,y1] (one realization)
% Recall y has cols [S,I,x]
data = vacc_sde( y0, tmax, k, choose_risk, wl, wh, d, num_realisations, sigma );


% Plot a realisation using tri_plot
my_title1=['Parameters: k = ',num2str(k), ...
    ', d = ',num2str(d),', sigma = ',num2str(sigma)];
my_title='';
total_plots = 10;
wcrit = d;
tri_plot(data,my_title,total_plots, wcrit)

