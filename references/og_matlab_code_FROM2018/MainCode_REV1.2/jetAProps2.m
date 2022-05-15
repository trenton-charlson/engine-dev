function [rho_c,viscK_c,cond_c, C_pc] = jetAProps2(Pc,Tc)
%jetAProps Solve for fluid properties of Jet-A given temperature and
%pressure - As of 04/18/18 no correlations are made to pressure but it will
%be included for future development as compressibility effects are
%considered

% Properties derived from 2 pt. linear correlation using DTIC A132106
% "Aviation Fuel Properties"  

%Solve for following properties:
% - Density
% - Viscosity
% - Conductivity
% - Specific heat

Tc_cel = Tc - 273.15; %Coolant temperature in celcius.

%% Density %%
% Temperature in Celcius
% rho_90C = 757 kg/m^3, rho_-40C = 848 kg/m^3, Pg. 22 
rhoSlope = -0.7000;
rho_c = rhoSlope*(Tc_cel + 40) + 848; %kg/m^3

%% Kinematic Viscosity %%
%Temperature in Kelvin
%From https://ws680.nist.gov/publication/get_pdf.cfm?pub_id=904848 -
%averaged between 2 referenced samples and corellated using linear
%interpolation
dat = load('jetA_visc.mat');
jetA_visc = dat.jetA_visc;
viscK_c = interp1(jetA_visc(:,1),jetA_visc(:,2),Tc,'linear').*(10^-6); %m^2/s

%% Specific Heat %%
%Temperature in Celcius
%Cp_40C = 2.04 kj/kgK, Cp_182C = 2.66 kj/kgK, Pg. 55
CpSlope = 0.0044;
C_pc = 1000.*(CpSlope.*(Tc_cel - 40) + 2.04); %J/kgK

%% Conductivity %%
%Temperature in Celcius
%Cond_20C = 0.1150 W/mK, Cond_100C = 0.1010 W/mK, Pg. 57
condSlope = -0.000175;
cond_c = condSlope.*(Tc_cel - 20) + 1.1150; %W/mK


end

