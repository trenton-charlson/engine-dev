function dP = pressureDrop2(CSA,perim,mDot_chan,rho_c,len,viscK_c)
%pressureDrop_passage - Calculate pressure drop in coolant passage
%% Initilization
epsilon = 0.005.*10^-3; %Surface roughness, guess based on our manufacuting capability, mm
d_hyd = 4.*CSA./perim; %Equivalent hydraulic diameter of coolant passage
w_c = (mDot_chan./(rho_c.*CSA)); %Coolant Velocity
Re_c = (w_c.*d_hyd)./viscK_c; %Reynolds number for coolant in passage (characteristic length chosen to be hydraulic diameter

%% Solve for Friction Factor
% Solved using Bellos, Nalbantis, Tsakiris friction estimation
%Valid for all flow regimes

a = 1/(1+(Re_c/2712)^8.4);
b = 1/(1+(Re_c/(150*d_hyd/epsilon)^1.8));

% Break up equation for troubleshooting
temp1 = (64/Re_c)^a;
temp2 = (.75*log(Re_c/5.37))^(2*b*(a-1));
temp3 = (.88*log(6.82*d_hyd/epsilon))^(2*(a-1)*(1-b));

fric = temp1*temp2*temp3;

%% Solve for pressure drop
% Pg. 93 Huzel and Huang 
dP = fric.*(len./d_hyd).*rho_c.*((w_c.^2)./(2*9.8)); %Solve for final dP across passage segment
    
end

