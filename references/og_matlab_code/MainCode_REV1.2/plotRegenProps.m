% T_cb = T_cb + dT_c;
%     T_wgvec(end+1) = T_wg;
%     T_cbvec(end+1) = T_cb;
%     visc_cvec(end+1) = viscK_cloc;
%     rho_cvec(end+1) = rho_cloc;
%     T_awvec(end+1) = T_aw;
%     T_wcvec(end+1) = T_wc;
%     qdot_vec(end+1) = qdot_ge;
%     h_gvec(end+1) = h_g;
%     vel_cvec(end+1) = vel_c;
%     P_cbvec(end+1) = P_cb;
%     P_cb = P_cb + dPc;
%     transitTime = transitTime + t_stay;

%% Plot T_cb, T_wg, T_wc %%
figure
subplot(2,1,1)
plot(Z.*1000,T_wgvec,'r')
hold on
plot(Z.*1000,T_cbvec,'b')
hold on
plot(Z.*1000,T_wcvec,'k')
title('T_{cb}, T_{wc}, T_{wg} vs. Axial Station')
ylabel('T, Kelvin')
legend('T_{wg}', 'T_{cb}', 'T_{wc}')

subplot(2,1,2)
plot(engineProps(:,2).*1000',engineProps(:,1).*1000','k')
hold on
plot(engineProps(:,2).*1000',zeros(1,numel(engineProps(:,2))),'-.r') %Centerline
xlabel('Axial Station Along Chamber, mm')
ylabel('Thrust Chamber Radius, mm')
xlim([0, 1100*engineProps(end,2)]);
%ylim([-5, (max(engineProps(:,1))*3500)]);

%% Plot vel_c, P_cb %%

figure
subplot(2,1,1)
plot(Z.*1000,vel_cvec,'k')
ylabel('Coolant Velocity, m/s')
hold on
yyaxis right
plot(Z.*1000,P_cbvec,'b')
ylabel('Pressure, Bar')
title('v_c & P_{cb} vs. Axial Station')
legend('v_c', 'p_{cb}', 'T_{wc}')

subplot(2,1,2)
plot(engineProps(:,2).*1000',engineProps(:,1).*1000','k')
hold on
plot(engineProps(:,2).*1000',zeros(1,numel(engineProps(:,2))),'-.r') %Centerline
xlabel('Axial Station Along Chamber, mm')
ylabel('Thrust Chamber Radius, mm')
xlim([0, 1100*engineProps(end,2)]);

%% Plot h_g, q_g %%

figure
subplot(2,1,1)
plot(Z.*1000,h_gvec,'k')
ylabel('Gas Side Convective Coeff., W/m^2K')
hold on
yyaxis right
plot(Z.*1000,qdot_vec,'b')
ylabel('Heat Flux, W/m^2')
title('h_g & q_g vs. Axial Station')
legend('h_g', 'q_g')

subplot(2,1,2)
plot(engineProps(:,2).*1000',engineProps(:,1).*1000','k')
hold on
plot(engineProps(:,2).*1000',zeros(1,numel(engineProps(:,2))),'-.r') %Centerline
xlabel('Axial Station Along Chamber, mm')
ylabel('Thrust Chamber Radius, mm')
xlim([0, 1100*engineProps(end,2)]);

%% Plot rho & visk %%

figure
subplot(2,1,1)
plot(Z.*1000,rho_cvec,'k')
ylabel('Coolant Density, kg/m^3')
hold on
yyaxis right
plot(Z.*1000,visc_cvec,'b')
ylabel('Coolant Viscosity, m^2/s')
title('\rho_c & \nu_c vs. Axial Station')
legend('\rho_c', '\nu_c')

subplot(2,1,2)
plot(engineProps(:,2).*1000',engineProps(:,1).*1000','k')
hold on
plot(engineProps(:,2).*1000',zeros(1,numel(engineProps(:,2))),'-.r') %Centerline
xlabel('Axial Station Along Chamber, mm')
ylabel('Thrust Chamber Radius, mm')
xlim([0, 1100*engineProps(end,2)]);








