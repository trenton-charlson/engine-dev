%% Generate Test Plots %%
%DATA FORMAT:
%R,Z,pip,aeat,mach,cf,ivac,isp,p,t,rho,h,u,mw,cp,gam,son,vis,cond,condfz,pran,pranfz,a_T,q_conv,T_aw
desProps = [10,9]; %Enter Values for desired property plots
propAxNames = {'Temperature, K','Pressure , Bar'};
plotNames = {'Nozzle Gas Temperature vs. Thrust Chamber Contour','Nozzle Static Pressure vs. Thrust Chamber Contour'};
contourShow = true; %Do or dont include chamber contour in each plot
for i = desProps
    figure
    if contourShow
        clear title xlabel ylabel
        plot(engineProps(:,2).*1000',engineProps(:,1).*1000','k') %Engine Contour
        xlabel('Axial Station Along Chamber, mm')
        title(plotNames{find(i == desProps)})
        hold on
        ylabel('Thrust Chamber Radius, mm')
        plot(engineProps(:,2).*1000',zeros(1,numel(engineProps(:,2))),'-.r') %Centerline
        xlim([0, 1100*engineProps(end,2)]);
        axis equal
        ylim([-5, (max(engineProps(:,1))*3500)]);
        axis manual
        hold on
        yyaxis right
        ylabel(propAxNames{find(i == desProps)})
        plot(engineProps(:,2).*1000',engineProps(:,i)','b')
        axis manual
    end
end