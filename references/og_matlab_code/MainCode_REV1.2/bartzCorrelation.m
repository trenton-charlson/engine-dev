function [h_g,q_conv,T_aw] = bartzCorrelation(engineProps,T_wg,R_t,P_c,C_star,R_tCurve,query)
%bartzCorrelation Calculates hot gas wall heat transfer coefficient for use
%in regen calcs. The Bartz method is employed.

%Query = points for which to run calc

%DATA FORMAT:
%R,Z,pip,aeat,mach,cf,ivac,isp,p,t,rho,h,u,mw,cp,gam,son,vis,cond,condfz,pran,pranfz
[numRows,~] = size(engineProps);


%% Calculate Parameters %%
%DUE TO BARTZ CORRELATION BEING DEVELOPED IN ENGLISH UNITS, VALUES
%CONVERTED TO ENGLISH UNITS IN EQUATIONS BELOW THEN RECONVERTED
g = 32.174; %Gravity, ft/s^2
T_wg = T_wg.*1.8; %Convert K to R
C_star = C_star.*3.28084; %Convert m/s to ft/s
R_tCurve = R_tCurve.*39.3701; %Convert m to in
R_t = R_t.*39.3701; %Convert m to in
P_c = P_c*14.5038; %Convert from Bar to PSI
T_c0 = engineProps(1,10).*1.8; % Convert K to R
gam = engineProps(query,16)'; %No conversion needed
gam0 = engineProps(1,16); %Chamber stagnation gamma
    %Old cp_ns - new estimate used as 2.2 KJ/Kg*K - correlation to come later
    % %NOTE: Effective specific heat used (instead of Frozen) 
    % cp_ns = engineProps(1,15).*0.23884; %Convert from Joules to BTU/lb*F
cp_ns = 2.2.*0.23884; %Convert from Kilo-Joules to BTU/lb*F
    % Originally used effective prandtl number from CEA, switched to gamma
    % correlation
    % praneff_ns = engineProps(1,22); %No conversion needed
praneff_ns = (4*gam0)./(9*gam0-5);
mach = engineProps(query,5)';
visc_ns = (engineProps(1,18)./1000).*(0.0672/12); %Convert to Poise then to lb/in-s
contourR = engineProps(query,1)'.*39.3701; %Convert to Inches
%Calculate Adiabatic Wall Temp
T_aw = T_c0.*(1+(praneff_ns.^(1/3)).*((gam-1)./2).*(mach.^2))./(1+((gam-1)./2).*(mach.^2));
%Sigma = sigA*sigB - Split into Terms
sigA = (0.5.*(T_wg./T_c0).*(1+((gam-1)./2).*mach.^2)+0.5).^(-0.68);
sigB = (1+((gam-1)./2).*(mach.^2)).^(-0.12);
sigma = sigA.*sigB;
%a_T = (aA*aB*aC*aD)*aE*sigma - Split into Terms
aA = 0.026./((2*R_t).^0.2);
aB = ((visc_ns.^0.2).*cp_ns)./(praneff_ns.^0.6);
aC = (P_c.*g./C_star).^0.8;
aD = ((2.*R_t)./R_tCurve).^0.1;
aE = ((R_t.^2)./(contourR.^2)).^0.9;
h_g = (aA.*aB.*aC.*aD).*aE.*sigma;
q_conv = h_g.*(T_aw-T_wg);

%% Re-Unit Convert %%
h_g = h_g.*144.*20428.175; %BTU/in^2-s-F -> BTU/ft^2-s-F -> W/m^2K
q_conv = q_conv.*1634246.235; %BTU/in^2-s -> W/m^2
T_aw = T_aw./1.8; %R -> K

% plot(engineProps(:,2)',h_g,'k');
% hold on
% yyaxis right
% plot(engineProps(:,2)',q_conv,'r');


end

