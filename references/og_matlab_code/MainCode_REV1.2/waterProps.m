function [rho_c,viscK_c,cond_c, C_pc] = waterProps(Pc,Tc)
%Solve water Properties for a given temperature
%Emperical Model Taken from "Equations of properties as a function of
%temperature for seven fluids"

% Density
rho_c = -3.0115*10^-6*Tc^3 + 9.6272*10^-4*Tc^2-.11052*Tc+1022.4; % Rho in kg/m^3

%Dynamic Viscosity
dynam = 3.8208*10^-2*(Tc-252.33)^-1; % Ns/m^2

%Kinematic Viscosity
viscK_c = dynam/rho_c ; %m^2/s

%Specific Heat
C_pc = 1000*(1.7850*10^-7*Tc^3-1.9149*10^-4*Tc^2+6.7953*10^-2*Tc-3.7559);  %J/kgK

%Conductivity
cond_c = 4.2365*10^-9*Tc^3-1.144*10^-5*Tc^2+7.1959*10^-3*Tc-.63262;% W/mk

end
