l:clear all
close all
nml=read_nml('params.nml','int_params.nml','modeselection.nml');
i=1;
energetics(i).name='zonal_heat_exchange';i=i+1;
energetics(i).name='oc_heat_exchange';i=i+1;
energetics(i).name='zonal_atm_LW_rec';i=i+1;
energetics(i).name='zonal_atm_LW_loss';i=i+1;
energetics(i).name='zonal_fric_ocean';i=i+1;
energetics(i).name='zonal_fric_internal';i=i+1;
energetics(i).name='zonal_pot2kin';i=i+1;
energetics(i).name='zonal_atm_SW_rec';i=i+1;
energetics(i).name='conv_zon2eddy_pot';i=i+1;
energetics(i).name='conv_zon2eddy_kin';i=i+1;
energetics(i).name='eddy_heat_exchange';i=i+1;
energetics(i).name='eddy_atm_LW_rec';i=i+1;
energetics(i).name='eddy_atm_LW_loss';i=i+1;
energetics(i).name='eddy_fric_ocean';i=i+1;
energetics(i).name='eddy_fric_internal';i=i+1; 
energetics(i).name='eddy_pot2kin';i=i+1;
energetics(i).name='ocean_bot_fric';i=i+1; 
energetics(i).name='ocean_wind_drag';i=i+1;
energetics(i).name='oc_LW_rad_rec';i=i+1;
energetics(i).name='oc_LW_rad_loss';i=i+1;
energetics(i).name='oc_SW_rad_rec';i=i+1;
energetics(i).name='zonalE_pot';i=i+1;
energetics(i).name='zonalE_kin';i=i+1;       
energetics(i).name='eddyE_pot';i=i+1;      
energetics(i).name='eddyE_kin';i=i+1;     
energetics(i).name='oceanE_thermal';i=i+1;    
energetics(i).name='oceanE_potential';i=i+1;  
energetics(i).name='oceanE_kinetic';

rawDATAmean=importdata('energetics_mean.dat');


for i=1:28
   % energetics(i).value=rawDATAmean(i,1);    
    energetics(i).mean=rawDATAmean(i,2);    
    energetics(i).xpower2=rawDATAmean(i,3);    
    energetics(i).xpower3=rawDATAmean(i,4);    
    energetics(i).xpower4=rawDATAmean(i,5);    
    energetics(i).cumulant3=rawDATAmean(i,6);    
    energetics(i).cumulant4=rawDATAmean(i,7);    
    energetics(i).variance=rawDATAmean(i,8);    
    energetics(i).skewness=rawDATAmean(i,9);    
    energetics(i).kurtosis=rawDATAmean(i,10);
end

% Function for displaying mean of a certain variable

meanof=@(str) energetics(cell2mat(arrayfun(@(x) any(strfind(x.name,str)),energetics,'UniformOutput',false))).mean;
varianceof=@(str) energetics(cell2mat(arrayfun(@(x) any(strfind(x.name,str)),energetics,'UniformOutput',false))).variance;
skewnessof=@(str) energetics(cell2mat(arrayfun(@(x) any(strfind(x.name,str)),energetics,'UniformOutput',false))).skewness;
kurtosisof=@(str) energetics(cell2mat(arrayfun(@(x) any(strfind(x.name,str)),energetics,'UniformOutput',false))).kurtosis;

