%% Regen Development Run %%
% This is the regen development code for engine analyis. It will use chamber
% properties generated by the main code to evaluate regen circuit seperate
% from primary engine analysis. Eventaually, it will be integrated in the
% main code deck.

%% Solution Architecture %%
% 0 = Inlet, 1 = Outlet
% 1. Solve for initial guess at circuit P0 & dP (P1-P0) using known P1 
%condition and assuming STP density
% 2. Solve for T1 using P0 guess, solve for actual P1
% 3. Guess new P0 based on previously calculated P1, repeat step 2 until
% results converge

%% Input Parameters - run dataCollectionScript.m first%%
fuelStiffness = .29; % 20% injector stiffness
P_ce_req = (P_c/(1-fuelStiffness))*bar;  %Required Coolant exit pressure, Pa
T_ci = 280.15; %Coolant Inlet Temperature, K
nstart = 102; %% Querey Point
%DATA FORMAT:
%R,Z,pip,aeat,mach,cf,ivac,isp,p,t,rho,h,u,mw,cp,gam,son,vis,cond,condfz,pran,pranfz,a_T,q_conv,T_aw
channelFile = 'channelDesign_YJ1S-03.xlsx';
channelProfile = channelProfileGen(channelFile);
channelProfile = channelProfile(1:nstart,:);
%channelProfile Layout
% R, Z, R_floor, R_ceil, numChannels, channelWidth, finWidth, channelDepth, 
%channelCSA, channelPerim, channelAR
engineProps=engineProps(1:nstart,:);
[numRows, numCols] = size(engineProps);
channelCSA = channelProfile(:,9);
numChannels = channelProfile(:,5);
channelPerim = channelProfile(:,10);
% mDot_chan = mDot_f./numChannels;
mDot_chan =2.8 ./numChannels;
Z = channelProfile(:,2);
R_floor = channelProfile(:,3);
R_ceil = channelProfile(:,4);

enginePerf = [R_t, P_c, C_star, R_tCurve]; %Variable Pass vector
%Will have to be reinvestigated
P_cb = 0; %bullshit initialization for bulk inlet pressure, doesn't yet affect properties
T_cb = T_ci; %Set coolant bulk temp to inlet temp
T_wgvec = [];
T_cbvec = [];
T_awvec = [];
T_wcvec = [];
qdot_vec = [];
h_gvec = [];
visc_cvec = [];
rho_cvec = [];
vel_cvec = [];
P_cbvec = [];
etavec =[];
hcvec =[];
transitTime = 0; %initialize clock for coolant transit
for i = numRows:-1:2 % Run heat transfer solution @ station i
    %Iterate heat transfer solution working from chamber exit to Injector
    [rho_cloc,viscK_cloc,cond_cloc, C_pcloc] = waterProps(P_cb,T_cb); %Solve for local props
    localProps = [rho_cloc,viscK_cloc,cond_cloc, C_pcloc];
    len = sqrt(((Z(i) - Z(i-1)).^2)+((R_floor(i) - R_floor(i-1)).^2));
    [dT_c, T_wc, T_wg, T_aw, h_g, qdot_ge, vel_c, t_stay,eta_tot,h_c] = HTran1D_Solve(T_cb,P_cb,engineProps,channelProfile,localProps,len,mDot_chan(i),enginePerf,i);
    dPc = pressureDrop2(channelCSA(i),channelPerim(i),mDot_chan(i),rho_cloc,len,viscK_cloc);
    T_cb = T_cb + dT_c;
    T_wgvec(end+1) = T_wg;
    T_cbvec(end+1) = T_cb;
    visc_cvec(end+1) = viscK_cloc;
    rho_cvec(end+1) = rho_cloc;
    T_awvec(end+1) = T_aw;
    T_wcvec(end+1) = T_wc;
    qdot_vec(end+1) = qdot_ge;
    h_gvec(end+1) = h_g;
    vel_cvec(end+1) = vel_c;
    P_cbvec(end+1) = P_cb;
    etavec(end+1)=eta_tot;
    hcvec(end+1)=h_c;
    P_cb = P_cb + dPc; % In Pascals
    transitTime = transitTime + t_stay;
end
%% Reverse them vectors to match coordinate system %%
T_cbvec(end+1) = T_cbvec(end) + (T_cbvec(end)-T_cbvec(end-1)); %apply dT_cb to end of coolant
T_cbvec = fliplr(T_cbvec); %reorient to match engine WCS
T_wgvec(end+1) = T_wgvec(end); %apply N-1 wall temp to end of iteration
T_wgvec = fliplr(T_wgvec);
visc_cvec(end+1) = visc_cvec(end); 
visc_cvec = fliplr(visc_cvec);
rho_cvec(end+1) = rho_cvec(end); 
rho_cvec = fliplr(rho_cvec);
T_awvec(end+1) = T_awvec(end); 
T_awvec = fliplr(T_awvec);
T_wcvec(end+1) = T_wcvec(end); 
T_wcvec = fliplr(T_wcvec);
qdot_vec(end+1) = qdot_vec(end); 
qdot_vec = fliplr(qdot_vec);
h_gvec(end+1) = h_gvec(end); 
h_gvec = fliplr(h_gvec);
vel_cvec(end+1) = vel_cvec(end); 
vel_cvec = fliplr(vel_cvec);
P_cbvec(end+1) = P_cbvec(end); 
P_cbvec = (P_cbvec + P_ce_req)/bar; %convert to Bar
etavec(end+1) = etavec(end);
etavec = fliplr(etavec);
hcvec(end+1)=hcvec(end);
hcvec = fliplr(hcvec);


% 
% plot(Z, T_wgvec,Z,T_wcvec)
% legend('T hot gas' ,'T coolant')
% xlabel('Axial Station Along Chamber, mm')
% ylabel('Temperature K')
% plot(Z, T_cbvec)
% title('Coolant Temp Vs Axial Station')
% xlabel('Axial Station Along Chamber, mm')
% ylabel('Coolant Temp K')
plot(Z, hcvec)
% title('Fin Efficiency Vs Axial Station')
% xlabel('Axial Station Along Chamber, mm')
% ylabel('Eff %')
%     