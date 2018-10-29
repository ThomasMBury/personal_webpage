% Script to run simulation of vaccine model drawing parameter values from
% triangular distributions.

close all
clear
tStart = tic;


% Specify input variables for vacc_sde_sensitivity
y0 = [0.01,1e-4,0.99];
tmax = 100;
k = 1000;
choose_risk = 1; 
wl = 1e-4;              
wh = 5e-4;  
d = 5e-4;
num_realisations = 200; 
sigma = 0.01;


% Run vacc_sde_sensitivity
data = vacc_sde_sensi( y0, tmax, k, choose_risk, wl, wh, d, num_realisations, sigma );

% Filter the data to 400 components in the desired interval using
% time_filter.m
tinit = 0;
tfin = tmax;
num_comps = 400;

data_filt = time_filter(data,tinit,tfin,num_comps);


% Save the filtered data to csv file for statistical analysis in R
csvwrite('sim_data/simdata_tri_sensi.txt',data_filt)


% Make plots using triplot.m
my_title1=['Parameters: k = ',num2str(k), ...
    ', d = ',num2str(d),', sigma = ',num2str(sigma)];
my_title='';
total_plots = 4;
w_crit = d;

tri_plot(data_filt,my_title,total_plots, w_crit)


% To save figure use print('Plots/plot_','-dpdf','-r0')

% End time counter
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));