
% Simulation of a critical trajectory near E_1 using vacc_sde
% Any time-series with epidemics are discarded - now look for indicators in
% I

clear
close all
tStart = tic;


% Input variables - see vacc_sde.m for definitions
y0 = [0.01,0,0.99];
tmax = 60;
k = 1000;
choose_risk = 3;
wl = 1e-4;              
wh = 7e-4;  
d = 5e-4;
num_realisations = 50; 
sigma = 0.01;

% Critical time when threshold crossed
if choose_risk == 3
    tcrit = tmax*(d-wl)/(wh-wl);
end

% Run vacc_sde
data = vacc_sde( y0, tmax, k, choose_risk, wl, wh, d, num_realisations, sigma );

%-------------
% Don't include time-series that have epidemics before CT, i.e. I>Ithresh
%--------------

Ithresh=1e-5;
tcomps=size(data,1);
epiless=[]; % vector of realisation numbers that don't have epidemics (I<Ithresh) in 0<t<tcrit

for i=1:num_realisations
    I=data(1:tcomps*(tcrit/tmax),3*i+1);
    check_vec = I > Ithresh;
    if sum(check_vec) == 0
        epiless=[epiless,i];
    end
end

% filtered data set
data_noepi=zeros(tcomps,2+3*length(epiless));
data_noepi(:,[1,2])=data(:,[1,2]);
for k=1:length(epiless)
    data_noepi(:,3*k : 3*k+2) = data(:,3*epiless(k):3*epiless(k)+2);
end
    



% Filter the data to 400 components in the desired interval using
% time_filter.m
tinit = 0;
tfin = tmax;
num_comps = 400;

data_filt = time_filter(data_noepi,tinit,tfin,num_comps);


% Save the filtered data to csv file for statistical analysis in R
csvwrite('sim_data/simdata_noepi2.txt',data_filt)


% Make plots using triplot.m
my_title1=['Parameters: k = ',num2str(k), ...
    ', d = ',num2str(d),', sigma = ',num2str(sigma)];
my_title='';
total_plots = 5;
w_crit = d;

tri_plot(data_filt,my_title,total_plots, w_crit)


% To save figure use print('Plots/plot_','-dpdf','-r0')

% End time counter
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));



