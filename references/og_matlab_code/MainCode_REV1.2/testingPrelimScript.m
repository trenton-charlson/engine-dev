%Monologue of ideas
% Outputs - points for countour
%         - Thermal analyis (Seperate Module)
%         - Comparative analysis
%         - Engine Properties (Guess at ISP, Thrust, Design Points, Mass Flow Est)
%         - Other Analysis of Engine

%% Constant Parameters %%
Ru = 8314.46; %Universal Gas Constant - J/kmolK
%% Input Parameters %%
%Uses LOX/Kerosene Data Deck
F_T = 2500; %Nominal Thrust Desired (N)
P_amb = 1.013529; %Design Exit Pressure (PSI) - For Nozzle
P_c = 14; %Design Chamber Pressure (Bar)
MR = 2.2; %Design Mixture Ratio O/F
P_rat = P_amb./P_c; %Design Pressure Ratio for Engine

%% Preliminary Combustion Chamber Sizing %%
%% ONLY FOR USE IN BASIC ANALYSIS %%
%Size Throat via determining Mass Flow Rate
k_c = 1.16; %GRAB FROM DATA DECK
T_c = 3209.91; %GRAB FROM DATA DECK
MW_c = 20.726; %GRAB FROM DATA DECK
k_t = 1.22; %GRAB FROM DATA DECK
P_t = 8.03; %GRAB FROM DATA DECK
T_t = 2440; %GRAB FROM DATA DECK
MW_t = 20.684; %GRAB FROM DATA DECK
V_e_ideal = sqrt(((2*k_c)./(k_c-1))*((Ru*T_c)/(MW_c))*(1-P_rat^((k_c-1)/k_c)));
mDot_total = F_T./V_e_ideal;
mDot_Ox = (mDot_total/(MR + 1))*MR;
mDot_F = mDot_total - mDot_Ox;
A_t = (mDot_total/(P_t*10^5))*sqrt((Ru*T_t)/(MW_t*k_t))
D_t = sqrt((4*A_t)/pi)
R_t = D_t/2
M_exit_ideal = sqrt((2/(k_t-1))*((1/P_rat)^((k_t-1)/k_t)-1))
isent1 = ((k_t+1)/2)^(1/(k_t-1))
isent2 = P_rat^(1/k_t)
isent3 = (k_t+1)/(k_t-1)
isent4 = P_rat^((k_t-1)/k_t)
A_e = A_t/(isent1*isent2*sqrt(isent3*isent4))
eR = A_e/A_t;

a = [1 0 0]
if any(a)
    b = 10
end






