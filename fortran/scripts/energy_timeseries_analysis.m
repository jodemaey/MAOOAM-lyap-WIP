clear all
close all
cd /net/aquana.mi/group1/u234069/MAOOAM-lyap-WIP/fortran/
% Prepare energetics
nml=read_nml('params.nml','int_params.nml','modeselection.nml');

% Read raw data of unformatted output file

T=nml.T_RUN/nml.TW;
fid=fopen('energetics_ts.dat');
rawData=fread(fid,T*29,'real*8');
T=length(rawData)/29;
fclose(fid);
fid=fopen('energetics_ts.dat');
rawData=reshape(fread(fid,T*29,'real*8'),29,T);
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
ocean_E_kinetic=rawData(i,:);i=i+1;
ocean_pot_tend=rawData(i,:);

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
    

set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextFontSize',10)

sample=[T-1000:T];
sample_time=sample*nml.F0^-1/3600/24/365;
figatm=figure;
subplot(2,2,1)
plot(sample_time,zonal_kinetic_energy_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(zonal_kinetic_energy_budget)./(1:T),'r--','Displayname','Mean')
plot(sample_time,zonal_fric_internal(sample),'displayname','zonal fric internal')
plot(sample_time,zonal_fric_ocean(sample),'displayname','zonal fric ocean')
plot(sample_time,zonal_pot2kin(sample),'displayname','zonal pot2kin')
plot(sample_time,-conv_zon2eddy_kin(sample),'displayname','-conv zon2eddy kin')
l=legend('show');set(l,'box','off');set(l,'color','none'); set(l,'location','eastoutside'); set(l,'orientation','vertical'); 
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on
title('zonal kinetic energy budget')
xlabel('Years')
ylabel('W/m^2')

subplot(2,2,2)
plot(sample_time,eddy_kinetic_energy_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(eddy_kinetic_energy_budget)./(1:T),'r--','Displayname','Mean')
plot(sample_time,eddy_fric_internal(sample),'displayname','eddy fric internal')
plot(sample_time,eddy_fric_ocean(sample),'displayname','eddy fric ocean')
plot(sample_time,eddy_pot2kin(sample),'displayname','eddy pot2kin')
plot(sample_time,conv_zon2eddy_kin(sample),'Displayname','conv zon2eddy kin')
l=legend('show');set(l,'box','off');set(l,'color','none'); set(l,'location','eastoutside'); set(l,'orientation','vertical'); 
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on
title('eddy kinetic energy budget')
xlabel('Years')
ylabel('W/m^2')

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

l=legend('show');set(l,'box','off');set(l,'color','none'); set(l,'location','eastoutside'); set(l,'orientation','vertical'); 
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on
title('zonal potential energy budget')
xlabel('Years')
ylabel('W/m^2')

subplot(2,2,4)
plot(sample_time,eddy_potential_energy_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(eddy_potential_energy_budget)./(1:T),'r--','Displayname','Mean')
plot(sample_time,eddy_atm_LW_rec(sample),'Displayname','eddy atm LW rec');
plot(sample_time,eddy_atm_LW_loss(sample),'Displayname','eddy atm LW loss');
plot(sample_time,eddy_heat_exchange(sample),'Displayname','eddy heat exchange');
plot(sample_time,-eddy_pot2kin(sample),'Displayname','-eddy pot2kin');
plot(sample_time,conv_zon2eddy_pot(sample),'Displayname','conv zon2eddy pot');

l=legend('show');set(l,'box','off');set(l,'color','none'); set(l,'location','eastoutside'); set(l,'orientation','vertical'); 
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on
title('eddy potential energy budget')
xlabel('Years')
ylabel('W/m^2')

set(gcf, 'Position', get(0, 'Screensize'));
saveas(figatm,'atm_energy_budget','fig')
export_fig('atm_energy_budget.pdf',figatm)
export_fig('atm_energy_budget.png',figatm)

figatmE=figure;
subplot(2,2,1)
plot((1:T)*nml.TW*nml.F0^-1/3600/24/365,zonal_E_kin,'k-')
title('zonal kinetic energy')
xlabel('Years')
ylabel('J/m^2')

subplot(2,2,2)
plot((1:T)*nml.TW*nml.F0^-1/3600/24/365,eddy_E_kin,'k-')
title('eddy kinetic energy')
xlabel('Years')
ylabel('J/m^2')

subplot(2,2,3)
plot((1:T)*nml.TW*nml.F0^-1/3600/24/365,zonal_E_pot,'k-')
title('zonal potential energy')
xlabel('Years')
ylabel('J/m^2')

subplot(2,2,4)
plot((1:T)*nml.TW*nml.F0^-1/3600/24/365,eddy_E_pot,'k-')
title('eddy potential energy')
xlabel('Years')
ylabel('J/m^2')


set(gcf, 'Position', get(0, 'Screensize'));
saveas(figatmE,'atm_energy','fig')
export_fig('atm_energy.pdf',figatmE)
export_fig('atm_energy.png',figatmE)


ocean_E_mech= ocean_E_kinetic+ocean_E_potential;

ocean_kinetic_budget=ocean_wind_drag+ocean_bot_fric-ocean_pot_tend;
ocean_potential_budget=ocean_pot_tend;
ocean_thermal_budget=oc_heat_exchange+oc_LW_loss+oc_LW_rec+oc_SW_rec;

figoc=figure;
subplot(2,3,1)
plot(sample_time,ocean_kinetic_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(ocean_kinetic_budget)./(1:T),'r--','Displayname','Mean')
plot(sample_time,ocean_bot_fric(sample),'displayname','ocean bot fric')
plot(sample_time,ocean_wind_drag(sample),'displayname','ocean wind drag')
plot(sample_time,-ocean_pot_tend(sample),'displayname','ocean potential tend')
l=legend('show');set(l,'box','off');set(l,'color','none'); set(l,'location','eastoutside'); set(l,'orientation','vertical'); 
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on
title('ocean kinetic energy budget')
xlabel('Years')
ylabel('W/m^2')

subplot(2,3,4)
plot((1:T)*nml.TW*nml.F0^-1/3600/24/365,ocean_E_kinetic,'k-')
title('ocean kinetic energy')
xlabel('Years')
ylabel('J/m^2')


subplot(2,3,2)
plot(sample_time,ocean_potential_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(ocean_potential_budget)./(1:T),'r--','Displayname','Mean')
l=legend('show');set(l,'box','off');set(l,'color','none'); set(l,'location','eastoutside'); set(l,'orientation','vertical'); 
set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on
title('ocean potential energy  budget')
xlabel('Years')
ylabel('W/m^2')

subplot(2,3,5)
plot((1:T)*nml.TW*nml.F0^-1/3600/24/365,ocean_E_potential,'k-')
title('ocean potential energy')
xlabel('Years')
ylabel('J/m^2')


subplot(2,3,3)
plot(sample_time,ocean_thermal_budget(sample),'k-','displayname','Budget')
hold on;
plot((1:T)*nml.F0^-1/3600/24/365,cumsum(ocean_thermal_budget)./(1:T),'r--','Displayname','Mean')
plot(sample_time,oc_heat_exchange(sample),'displayname','oc heat exchange')
plot(sample_time,oc_LW_loss(sample),'displayname','oc LW loss')
plot(sample_time,oc_LW_rec(sample),'displayname','oc LW rec')
plot(sample_time,oc_SW_rec(sample),'Displayname','oc SW rec')
l=legend('show');set(l,'box','off');set(l,'color','none'); set(l,'location','eastoutside'); set(l,'orientation','vertical');

set(gca,'xlim', [sample_time(1) sample_time(end)])
grid on
title('ocean thermal budget')
xlabel('Years')
ylabel('W/m^2')

subplot(2,3,6)
plot((1:T)*nml.TW*nml.F0^-1/3600/24/365,ocean_E_thermal,'k-')
title('ocean thermal energy')
xlabel('Years')
ylabel('J/m^2')

set(gcf, 'Position', get(0, 'Screensize'));

saveas(figoc,'oc_energy','fig')
export_fig('oc_energy.pdf',figoc)
export_fig('oc_energy.png',figoc)
