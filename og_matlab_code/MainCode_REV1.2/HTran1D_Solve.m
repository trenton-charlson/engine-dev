function [dT_c, T_wc, T_wg, T_aw, h_g, qdot_ge, vel_c, t_stay,eta_tot,h_c] = HTran1D_Solve(T_cb,P_cb,engineProps,channelProfile,localProps,len,mDot_chan,enginePerf,i)
%1DHT_Solve - 1D Heat transfer analysis code for YJ-1S Development. Solves
%for thermal equilibrium @ station I

%channelProfile Layout
% R, Z, R_floor, R_ceil, numChannels, channelWidth, finWidth, channelDepth, 
%channelCSA, channelPerim, channelAR

%engineProps Layout:
%R,Z,pip,aeat,mach,cf,ivac,isp,p,t,rho,h,u,mw,cp,gam,son,vis,cond,condfz,pran,pranfz,a_T,q_conv,T_aw,T_wg

%i = engine station point. T_cbi = coolant inlet temperature into station.
%P_cbi = coolant inlet pressure into station

%% Material Parameters & Constants %%
cond_w = 167; %W/m-k, wall conductivity for 6061

%% Input Common Variables %%
%Channel Properties
R = channelProfile(i,1);
t_fin = channelProfile(i,7); %fin thickness
N = channelProfile(i,5); %number of channels
Z = channelProfile(i,2);
R_floor = channelProfile(i,3);
R_ceil = channelProfile(i,4);
a = channelProfile(i,6); %channel width
L = channelProfile(i,8); %channel depth
A_cc = channelProfile(i,9); %coolant channel CSA
perim_hyd = channelProfile(i,10);%channel hydraulic perimeter

%Coolant Properties
rho_c = localProps(1); %Coolant Density kg/m^3
viscK_c = localProps(2); %Coolant viscosity m^2/s
cond_c = localProps(3); %Coolant conductivity W/mK
C_pc = localProps(4); %Coolant specific heat J/kgK

%Engine Geometry, Etc.
R_t = enginePerf(1);
P_c = enginePerf(2);
C_star = enginePerf(3);
R_tCurve = enginePerf(4);

%% Solver Rev 2 %%
%len = station fin width

converged = false; %convergence flag to exit solver loop
iter = 1; %iteration count

%Geometric Parameters (AT STATION, NOT SUM)
D_hyd = (4*A_cc)/perim_hyd; %Channel Hydraulic diameter
L_cor = L + t_fin/2; %Corrected channel depth for fin efficiency calcs
perim_conf = 2*len + 2*t_fin; %contact perimeter for fin efficiency calcs
A_conf = len*t_fin; %contact area for fin efficiency calcs
A_fin = 2*L_cor*len; %area of single fin in analysis segment, corrected for adiabatic tip
A_finTot = N*A_fin; %total fin area 
A_totc = A_finTot + N*a*len; %total coolant side heat transfer area 
A_totg = 2*pi*R*len; %total hot gas side wall area for heat transfer

%Coolant Side Parameters
vel_c = mDot_chan/(A_cc*rho_c); %coolant velocity in channel (m/s)
t_stay = len/vel_c; %coolant stay time in channel segment
Re_c = (vel_c.*D_hyd)./viscK_c; %Coolant Reynolds number
Pr_c = viscK_c.*rho_c.*C_pc./cond_c; %Coolant Prandtl Number

%Calculate Gas Side Estimate
T_wg1 = 800; %initial gas side temp estimate
fudge = 1;

while ~converged % Loop till steady state solution
    %Solve for gas side heat flux estimate
    [h_g,qdot_ge,T_aw] = bartzCorrelation(engineProps,T_wg1,R_t,P_c,C_star,R_tCurve,i);
    q_ge = qdot_ge*A_totg*fudge;
    T_wc1 = T_wg1 - (qdot_ge*(R_floor-R)/cond_w);
    h_c = 0.021*(Re_c^0.8)*(Pr_c^0.4)*(0.64+0.36*(T_cb/T_wc1))*(cond_c/D_hyd); %coolant side
    m = sqrt((h_c*perim_conf)/(cond_w*A_conf)); %fin efficiency parameter
    %Ref. Fundamentals of heat and mass transfer, 6e, pp. 152, square fin
    %coefficient for adiabatic tip condition
    eta_fin = tanh(m*L_cor)/(m*L_cor); %fin efficiency
    eta_tot = 1-(((N*A_fin)/A_totc)*(1-eta_fin)); %overall coolant side heat transfer efficiency
    q_ce = h_c*A_totc*eta_tot*(T_wc1-T_cb); %Coolant side heat transfer estimate
    
    %Begin convergence looping
    T_wg2 = T_wg1; %set up storage variable for previous gas value
    err = abs((q_ce-q_ge)/q_ge);
    if err < 0.005
        converged = true;
    elseif iter > 100
        converged = true;
        %Print "error" message:
        'Solution took 100 Iterations - moving on'
        i
    else
       if q_ge > q_ce
           T_wg1 = T_wg1+(100*err); %make next guess at T_wg
       elseif q_ce>q_ge
           T_wg1 = T_wg1-(100*err);
       end
    end
    iter = iter + 1;
    
end

dT_c = ((q_ce+q_ge)/2)/(C_pc*mDot_chan*N);

T_wg = T_wg1;
T_wc = T_wc1;


%% OLD CODE %%
%% Solver 1 %%
% A_wg = ((channelFloorArcLn(i) + channelFloorArcLn(i-1))./2).*len; %m^2
% %Solve for hot gas area using average channel floor arc across inlet and outlet of
% %station
% d_hyd = 4.*channelCSA(i)./perim(i); %Equivalent hydraulic diameter of coolant passage, m
% w_c = (mDot_chan(i)./(rho_ci.*channelCSA(i))); %Coolant velocity
% q_convi = engineProps(i,24); %Guess at Gas side convective heat flux (from prev iteration)
% T_wgi = engineProps(i,26); %Guess at T_wg (from prev iteration)
% Re_c = (w_c.*d_hyd)./viscK_ci; %Coolant Reynolds number
% Pr_c = viscK_ci.*rho_ci.*C_pci./cond_ci; %Prandtl Number
% dT_c = (q_convi.*A_wg)./(mDot_chan(i).*localProps(4)); %Calculated coolant temperature gain
% T_wci = T_wgi - (q_convi.*(R_floor(i)-R(i))./cond_w); %GasHT Coolant side wall temp
% %Calculate T_wc based on wall conductivity at temperature
% Nu_c = 0.021.*(Re_c.^0.8).*(Pr_c.^0.4).*(0.64 + 0.36.*(T_cbi./T_wci)); %Calculated coolant Nusselt Number
% a_c = Nu_c.*(cond_ci./d_hyd); %Calculated coolant heat transfer coefficient
% T_wc_new = (q_convi./a_c) + T_cbi; %Calculate new guess at T_wc based on coolant side HT Coeff
% T_wg_new = (q_convi.*(R_floor(i)-R(i))./cond_w) + T_wc_new; %Calculate new guess at T_wg based on T_wc guess





% TwgErr = abs((T_wgi - T_wg_new)./T_wgi);
% if TwgErr < 0.05
%     converged = true;
% else
%     converged = false;
% end

end

