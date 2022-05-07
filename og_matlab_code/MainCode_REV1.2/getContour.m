function [engineContour,chBarrel,nozzleContour, R_tCurve] = getContour(R_t,LStar,conRat,conAng,radRat,expRat,plot1)
%Generate Nozzle Contour given the following input parameters:
% R_t - Throat Radius (m)
% LStar - L* Value (derived from Propellants)
% conRat - Contraction Ratio for chamber
% conAng - Convergence Angle for chamber
% radRat - Radius Ratio R/R_cyl for chamber convergence section lead in radius,
% where R_cyl is the cylindrical chamber radius
%% Size Chamber Parameters %%
exitHalfAngle = 15; %Nozzle Exit half angle, degrees
throatLeadInRadius = 1.5*R_t; %Throat lead in radius after contraction angle
throatLeadOutRadius = 0.382*R_t; %Nozzle lead in radius
R_tCurve = (throatLeadInRadius+throatLeadOutRadius)/2; %Throat radius of curvature (for Bartz)
throatArea = pi*R_t^2; %Area of the throat (duh)
chamberArea = conRat * throatArea;
chamberRad = sqrt(conRat)*R_t;
%% Size Nozzle Critical Points %%
% (z – h)^2 + (r – k)^2 = r^2
%Critical Pts:
%     - 0: Chamber Start, Injector Plane
%     - 1: Converging radius start
%     - 2: Converging radius end, angular section start
%     - 3: Angular section end, nozzle lead in start
%     - 4: Throat - Start of nozzle lead in (minor contour)
%     - 5: Start of major nozzle contour
%     - 6: End of Nozzle
%Size Angular Section Lead in
throatLeadRDist = throatLeadInRadius-(cosd(conAng)*throatLeadInRadius)+R_t; %R distance occupied by throat lead in
convLdRDelta = chamberRad - throatLeadRDist; %Delta R for converging lead in
conLdRadMax = convLdRDelta/(1-cosd(conAng)); %Max lead in rad for meet @conAng
conLeadInRadius = radRat*conLdRadMax;
a1r = chamberRad - conLeadInRadius; %Radial center point for lead out arc
r0 = chamberRad;
%Chamber lead out start
z1 = 0;
r1 = chamberRad;
%Chamber lead out end, Angular contraction start
z2 = z1+conLeadInRadius*sind(conAng); %axial coordinate of end of lead in rad
r2 = chamberRad - conLeadInRadius*(1-cosd(conAng)); %Set radial height of ending coord
a2r = throatLeadInRadius + R_t; %Radial center point for lead in
%Angular contraction end,throat lead in start
r3 = throatLeadRDist;
z3 = (r3-r2)/(-tand(conAng)) + z2;
a3r = throatLeadOutRadius + R_t;
%Throat
z4 = z3+throatLeadInRadius*sind(conAng);
r4 = R_t;
%Throat lead out start
z5 = z4+(throatLeadOutRadius*sind(exitHalfAngle));
r5 = throatLeadOutRadius - sqrt(throatLeadOutRadius^2-(throatLeadOutRadius*sind(exitHalfAngle))^2) + R_t;
%Nozzle Sizer
%Currently conic only - mabye parabolic in the future.
r6 = sqrt(expRat)*R_t;
z6 = z5+(r6-r5)/(tand(exitHalfAngle));

%% Generate Nozzle Contour Vector %%
nozzleLength = z6; %Total Length of Nozzle
numPTS = 100;
nozzleZ = linspace(0,nozzleLength,numPTS);
nozzleR = [];
for i = 1:numel(nozzleZ)
    if z1 <= nozzleZ(i) && nozzleZ(i) < z2
        %Chamber Lead Out Rad
        nozzleR(i) = (chamberRad - conLeadInRadius)+sqrt(conLeadInRadius^2-nozzleZ(i)^2);
    elseif z2 <= nozzleZ(i) && nozzleZ(i) < z3
        nozzleR(i) = r2-(nozzleZ(i)-z2)*tand(conAng);
    elseif z3 <= nozzleZ(i) && nozzleZ(i) < z4
        nozzleR(i) = R_t + throatLeadInRadius - sqrt(throatLeadInRadius^2-(z4-nozzleZ(i))^2);
    elseif z4 <= nozzleZ(i) && nozzleZ(i) < z5
        nozzleR(i) = R_t + throatLeadOutRadius - sqrt(throatLeadOutRadius^2-(nozzleZ(i)-z4)^2);        
    elseif z5 <= nozzleZ(i) && nozzleZ(i) <= z6
        nozzleR(i) = r5+(nozzleZ(i)-z5)*tand(exitHalfAngle);
    else
        nozzleR(i) = 0;
    end    
end
%plot(nozzleZ,nozzleR,'k')
% hold on

%% Size Combustion Chamber %%
chamberVolume = throatArea*LStar;
%Solve for convergent volume:
convergingVolume = 0;
[nozzleZ, sortInd] = sort([nozzleZ,z4]);
nozzleR = [nozzleR,r4];
nozzleR = nozzleR(sortInd);
[~,throatInd] = max(sortInd);
for i = 2:throatInd
    deltaVol = pi * (nozzleZ(i)-nozzleZ(i-1)) * mean([nozzleR(i),nozzleR(i-1)])^2;
    convergingVolume = convergingVolume + deltaVol;
end
cylVol = chamberVolume-convergingVolume;
cylLength = cylVol/(pi*(chamberRad^2))
chBarrelZ = linspace(0,cylLength,numPTS);
chBarrelR = ones(1,numel(chBarrelZ)).*chamberRad;

%% Shift Points and Generate Outputs
zCrit = [0 ([z1,z2,z3,z4,z5,z6] + cylLength)]; %Critical Z points per defs 1-6 above
rCrit = [r0,r1,r2,r3,r4,r5,r6];
chBarrel = [chBarrelR;chBarrelZ];%Points for chamber barrel
nozzleZ = nozzleZ + cylLength;
nozzleContour = [nozzleR;nozzleZ]; %Points for nozzle section
engineContour = [chBarrel,[nozzleR;(nozzleZ+cylLength)]];
%plot(engineContour(2,:),engineContour(1,:)) 


%% Generate Plots %%
% %Adapt Coordinate System
% nozzleRVector = [r45Vector,r56Vector];
% nozzleZVector = [z45Vector,z56Vector] + chBarrelZ(end); %Offset by Combustion Chamber
% zConvVector = zConvVector + chBarrelZ(end);
% rTotal = [chBarrelR(1:end-1),rConvVector(1:end-1),nozzleRVector];
% zTotal = [chBarrelZ(1:end-1),zConvVector(1:end-1),nozzleZVector];
% chamberContour = [zTotal;rTotal];
% chBarrel = [chBarrelR;chBarrelZ];
% convContour = [rConvVector(2:end);zConvVector(2:end)];
% divContour = [nozzleRVector(2:end);nozzleZVector(2:end)];
% 
% if plot1
% plot(zTotal,rTotal,'k')
% axis equal
% xlim([0, 1.1*nozzleZVector(end)]); %m
% ylim([-0.05*chamberRad, 1.1*chamberRad]); %m
% hold on
% plot(zTotal,zeros(1,numel(zTotal)),'-.r')
% title('Chamber Contour')
% xlabel('Z, mm')
% ylabel('R, mm')
% end

%% Plot Errything (DEBUGGING)%%
% plot(z12Vector, r12Vector)
% axis equal
% hold on
% plot(z34Vector, r34Vector)
% hold on
% plot(z23Vector, r23Vector)
% hold on 
% plot(z45Vector, r45Vector)
% hold on 
% plot(z56Vector, r56Vector)
% plot(zConvVector,rConvVector)
%slop = (r12Vector(end)-r12Vector(end-1))/(z12Vector(end)-z12Vector(end-1))


end

