clear all
close all
% Prepare energetics
nml=read_nml('params.nml','int_params.nml','modeselection.nml');

% Read raw data of unformatted output file

T=nml.T_RUN/nml.TW;
fid=fopen('energetics_ts.dat');
rawData=fread(fid,T*28,'real*8');
T=length(rawData)/28;
fclose(fid);
fid=fopen('energetics_ts.dat');
rawData=reshape(fread(fid,T*28,'real*8'),28,T);
fclose(fid);

i=1;
zonal_heat_exchange=rawData(i,:);i=i+1;
oc_heat_exchange=rawData(i,:);i=i+1;
zonal_atm_LW_rec=rawData(i,:);i=i+1;
zonal_atm_LW_loss=rawData(i,:);i=i+1;
zonal_fric_ocean=rawData(i,:);i=i+1;
zonal_fric_internal=rawData(i,:);i=i+1;
zonal_pot2kin=rawData(i,:);i=i+1;
zonal_atm_SW_rec=rawData(i,:);i=i+1;
conv_zon2eddy_pot=rawData(i,:);i=i+1;
conv_zon2eddy_kin=rawData(i,:);i=i+1;
eddy_heat_exchange=rawData(i,:);i=i+1;
eddy_atm_LW_rec=rawData(i,:);i=i+1;
eddy_atm_LW_loss=rawData(i,:);i=i+1;
eddy_fric_ocean=rawData(i,:);i=i+1;
eddy_fric_internal=rawData(i,:);i=i+1; 
eddy_pot2kin=rawData(i,:);i=i+1;
ocean_bot_fric=rawData(i,:);i=i+1; 
ocean_wind_drag=rawData(i,:);i=i+1;
oc_LW_rec=rawData(i,:);i=i+1;
oc_LW_loss=rawData(i,:);i=i+1;
oc_SW_rec=rawData(i,:);i=i+1;
zonal_E_pot=rawData(i,:);i=i+1;
zonal_E_kin=rawData(i,:);i=i+1;       
eddy_E_pot=rawData(i,:);i=i+1;      
eddy_E_kin=rawData(i,:);i=i+1;     
ocean_E_thermal=rawData(i,:);i=i+1;    
ocean_E_potential=rawData(i,:);i=i+1;  
ocean_E_kinetic=rawData(i,:);


%
% Derived quantities
%

atm_E_total=zonal_E_kin+zonal_E_pot+eddy_E_kin+eddy_E_pot;
zonal_kinetic_energy_budget= ...
    zonal_fric_internal+...
    zonal_fric_ocean + ...
    zonal_pot2kin+...
    -conv_zon2eddy_kin;
eddy_kinetic_energy_budget= ...
    eddy_fric_internal+...
    eddy_fric_ocean + ...
    eddy_pot2kin+...
    +conv_zon2eddy_kin;
total_kinetic_energy_budget= ...
    zonal_fric_internal+...
    zonal_fric_ocean + ...
    eddy_fric_internal+...
    eddy_fric_ocean+ ...
    zonal_pot2kin+ ...
    eddy_pot2kin;

zonal_potential_energy_budget= ...
    zonal_atm_LW_rec+...
    zonal_atm_LW_loss+...
    zonal_atm_SW_rec+...
    zonal_heat_exchange + ...
    -zonal_pot2kin+...
    -conv_zon2eddy_pot;
eddy_potential_energy_budget= ...
    eddy_atm_LW_rec+...
    eddy_atm_LW_loss+...
    eddy_heat_exchange + ...
    -eddy_pot2kin...
    +conv_zon2eddy_pot;
total_potential_energy_budget= ...
    zonal_atm_LW_rec+...
    zonal_atm_LW_loss+...
    zonal_atm_SW_rec+...
    zonal_heat_exchange + ...
    -zonal_pot2kin+...
    eddy_atm_LW_rec+...
    eddy_atm_LW_loss+...
    eddy_heat_exchange + ...
    -eddy_pot2kin;
    
    
% % Function for obtaining timeseries and so on
% cut=1; % This allows to chop off some initial points in each time series
% cutoff=@(x) x(cut:end); T=T-cut+1;
% tsof=@(str) cutoff(energetics(cell2mat(arrayfun(@(x) any(strfind(x.name,str)),energetics,'UniformOutput',false))==1).data);
% meanofts=@(x) mean(tsof(x));


sample=[T-1000:T];
sample_time=sample*nml.F0^-1/3600/24/365;
figure;
subplot(2,2,1)
plot(sample_time,zonal_kinetic_energy_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(zonal_kinetic_energy_budget)./(1:T),'r--','Displayname','Mean')
plot(sample_time,zonal_fric_internal(sample),'displayname','zonal fric internal')
plot(sample_time,zonal_fric_ocean(sample),'displayname','zonal fric ocean')
plot(sample_time,zonal_pot2kin(sample),'displayname','zonal pot2kin')
plot(sample_time,-conv_zon2eddy_kin(sample),'displayname','-conv zon2eddy kin')
legend('show')
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on
title('zonal_kinetic_energy_budget')

subplot(2,2,2)
plot(sample_time,eddy_kinetic_energy_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(eddy_kinetic_energy_budget)./(1:T),'r--','Displayname','Mean')
plot(sample_time,eddy_fric_internal(sample),'displayname','eddy fric internal')
plot(sample_time,eddy_fric_ocean(sample),'displayname','eddy fric ocean')
plot(sample_time,eddy_pot2kin(sample),'displayname','eddy pot2kin')
plot(sample_time,conv_zon2eddy_kin(sample),'Displayname','conv zon2eddy kin')
legend('show')
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on

title('eddy_kinetic_energy_budget')

subplot(2,2,3)
plot(sample_time,zonal_potential_energy_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(zonal_potential_energy_budget)./(1:T),'r--','Displayname','Mean')
plot(sample_time,zonal_atm_LW_rec(sample),'Displayname','zonal atm LW rec');
plot(sample_time,zonal_atm_LW_loss(sample),'Displayname','zonal atm LW loss');
plot(sample_time,zonal_atm_SW_rec(sample),'Displayname','zonal atm SW rec');
plot(sample_time,zonal_heat_exchange(sample),'Displayname','zonal heat exchange');
plot(sample_time,-zonal_pot2kin(sample),'Displayname','-zonal pot2kin');
plot(sample_time,-conv_zon2eddy_pot(sample),'Displayname','-conv zon2eddy pot');

legend('show')
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on
title('zonal_potential_energy_budget')

subplot(2,2,4)
plot(sample_time,eddy_potential_energy_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(eddy_potential_energy_budget)./(1:T),'r--','Displayname','Mean')
plot(sample_time,eddy_atm_LW_rec(sample),'Displayname','eddy atm LW rec');
plot(sample_time,eddy_atm_LW_loss(sample),'Displayname','eddy atm LW loss');
plot(sample_time,eddy_heat_exchange(sample),'Displayname','eddy heat exchange');
plot(sample_time,-eddy_pot2kin(sample),'Displayname','-eddy pot2kin');
plot(sample_time,conv_zon2eddy_pot(sample),'Displayname','conv zon2eddy pot');

legend('show')
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on
title('eddy_potential_energy_budget')


sample=1:T;
sample_time=sample*nml.F0^-1/3600/24/365;
figure;
subplot(2,2,1)
plot(sample_time,zonal_kinetic_energy_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(zonal_kinetic_energy_budget)./(1:T),'r--','Displayname','Mean')
legend('show')
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on
subplot(2,2,2)
plot(sample_time,eddy_kinetic_energy_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(eddy_kinetic_energy_budget)./(1:T),'r--','Displayname','Mean')
legend('show')
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on

subplot(2,2,3)
plot(sample_time,zonal_potential_energy_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(zonal_potential_energy_budget)./(1:T),'r--','Displayname','Mean')

legend('show')
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on

subplot(2,2,4)
plot(sample_time,eddy_potential_energy_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(eddy_potential_energy_budget)./(1:T),'r--','Displayname','Mean')

legend('show')
grid on


sample=1:T;
sample_time=sample*nml.F0^-1/3600/24/365;
figure;
subplot(2,2,1)
%plot(sample_time,zonal_kinetic_energy_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(zonal_kinetic_energy_budget)./(1:T),'r--','Displayname','Mean')
legend('show')
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on
subplot(2,2,2)
plot(sample_time,zonal_kinetic_energy_budget(sample)+eddy_kinetic_energy_budget(sample)+zonal_potential_energy_budget(sample)+eddy_potential_energy_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(zonal_kinetic_energy_budget(sample)+eddy_kinetic_energy_budget(sample)+zonal_potential_energy_budget(sample)+eddy_potential_energy_budget(sample))./(1:T),'r--','Displayname','Mean')
legend('show')
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on

subplot(2,2,3)
plot(sample_time,zonal_kinetic_energy_budget(sample)+eddy_kinetic_energy_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(zonal_kinetic_energy_budget(sample)+eddy_kinetic_energy_budget(sample))./(1:T),'r--','Displayname','Mean')

legend('show')
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on

subplot(2,2,4)
plot(sample_time,zonal_potential_energy_budget(sample)+eddy_potential_energy_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(zonal_potential_energy_budget(sample)+eddy_potential_energy_budget(sample))./(1:T),'r--','Displayname','Mean')

legend('show')
grid on
