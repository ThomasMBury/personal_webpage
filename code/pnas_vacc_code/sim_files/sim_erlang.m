% Simulation of a critical trajectory near E_1 using vacc_sde_erlang

clear
close all
tStart = tic;


% Input variables - see vacc_sde.m for definitions
y0 = ones(22,1)*1e-5; y0(1)=0.01; y0(22)=0.99;
tmax = 100;
k = 1000;
choose_risk = 1;
wl = 1e-4;              
wh = 5e-4;  
d = 5e-4;
num_realisations = 200; 
sigma = 0.01;

% Critical time when threshold crossed
if choose_risk == 3
    tcrit = tmax*(d-wl)/(wh-wl);
end

% Run vacc_sde
data = vacc_sde_erlang_seasonal( y0, tmax, k, choose_risk, wl, wh, d, num_realisations, sigma );

% Filter the data to 400 components in the desired interval using
% time_filter.m
tinit = 0;
tfin = tmax;   
num_comps = 400;

data_filt = time_filter(data,tinit,tfin,num_comps);


% Save the filtered data to csv file for statistical analysis in R
csvwrite('sim_data/simdata_erlangseasonal_highkap_tri.txt',data_filt)
% save('sim_data/simdata_temp.mat','data_filt')


% Make a figure
t=data_filt(:,1);
figure()
plot_num=1;
yyaxis left
plot(t,data_filt(:,22*plot_num+2),'b')
axis([-inf,inf,0.5,1])
set(gca,'FontSize',12,'YColor','b')

yyaxis right
plot(t,sum(data_filt(:,22*plot_num-18 : 22*plot_num+1),2),'r')
set(gca,'YColor','r')
axis([-inf,inf,-inf,inf])



% To save figure use print('Plots/plot_','-dpdf','-r0')

% End time counter
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));



