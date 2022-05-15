%% Data Import INFO %%
% Goal of importing data from excel deck. Data to be imported is a user
% generated CEA output deck and can contain data across Contraction ratios,
% seperated by sheets, and data for any number of chamber pressures. Data
% for a given chamber pressure seperates each run within one contraction
% ratio group. Each chamber pressure datum includes combustion calculations
% at the chamber, combustion end, throat, up to 13 user defined Pi/Pe
% ratios, and up to 13 user defined supersonic expansion area ratios
% (Ae/At) 

%% Data Import HEADER %%
% The data from each each spreadsheet shall have the following header and
% column format:
% cRat,pip,aeat,mach,cf,ivac,isp,p,t,rho,h,u,mw,cp,gam,son,vis,cond,condfz,pran,pranfz

%% Constants %%
g0 = 9.81; %Gravity, m/s^2
bar = 100000; %1 Bar in Pa
Ru = 8314.46; %Universal Gas Constant - J/kmolK

%% Specify Parameters %%
thrust_nom = 3500; %Nominal design thrust in Newtons
P_c = 18; %Chamber Pressure in Bar
P_amb = 1.01325; %Ambient Pressure in Bar
conRat = 6; %Contraction Ratio for Engine
P_rat1 = P_amb/P_c; %Pressure ratio of overall expansion
P_rat2 = P_c/P_amb;

%% Propellants: Kerosene and Liquid Oxygen
MR = 1.8; %Mixture Ratio: MUST MATCH MR USED IN DATA DECK GENERATION
LStar = 1.05; %Characteristic length, meters

%% Data Retrieval and Interpolation %%
fileName = 'dataDeck_Rev1_MR1.8_cR3-8.xlsx';
conRat_array = [3, 4, 5, 6, 7, 8]; %Contraction Ratios from Data Deck
sheetCount = numel(conRat_array); %Number of sheets to pull data from
[~,sheet_name]=xlsfinfo(fileName); %pull datasheet names in
pStations = [6,8,10,12,14,16,18,20]; %Pressure Stations for Data Deck
pRatStations = [1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.5,1.6,1.7,1.8,1.9]; % Pressure ratio stations for data deck (Converging)
areaRatStations = [1.1,1.2,1.3,1.5,1.7,1.9,2.1,2.3,2.6,3,3.4,3.8,4.2]; % Area Ratio stations for Data Deck (Diverging)
import = true; %Toggle whether or not dataMain must be regenerated
if sheetCount ~= numel(sheet_name)
    error = 'ERROR, number of sheets does not match contraction ratio count'
    return;
end
if import
    dataMain = {};
    for k=1:numel(sheet_name) % Establish data matrix for given contraction ratios;
      dataMain{k}=xlsread(fileName,sheet_name{k});
    end
end
%Get Data Deck for Contour Analysis
CEAData_int = getDataDeck(dataMain,P_c,conRat,conRat_array,pStations); %Get CEA Data Deck via interpolation
[mainProperties,nozzleProps_working] = getNozzleData(P_rat2,conRat,CEAData_int,pRatStations,areaRatStations); %Interpolate/Derive Nozzle Properties based on Area Ratio

%% Size & Calculate Engine Parameters %%
% DATA LAYOUT: 
%pip,aeat,mach,cf,ivac,isp,p,t,rho,h,u,mw,cp,gam,son,vis,cond,condfz,pran,pranfz
divAngle = 15; %Divergence half angle, degrees
nozCorrFactor = (1+cosd(divAngle))/2; %Correction factor for Nozzle, using conic appx
V_exit = mainProperties(4,3)*mainProperties(4,15); %Sonic Vel * Mach @ Exit
Cf_exit = mainProperties(4,4); %Exit thrust coeficient
mDot_tot = thrust_nom/(nozCorrFactor*V_exit); %Calculate total mass flow based on V_exit(kg/s)
mDot_o = (mDot_tot/(MR + 1))*MR; %Oxidizer mass flow
mDot_f = mDot_tot - mDot_o; %fuel mass flow

%% Size Throat, Chamber, and Nozzle %%
T_c = mainProperties(1,8);
A_t = (mDot_tot/(mainProperties(3,7)*bar))*sqrt((Ru*mainProperties(3,8))/(mainProperties(3,12)*mainProperties(3,14))); %Calculate Throat Area (m^2)
C_star = (P_c.*bar.*A_t./mDot_tot)*0.975; %Calculate C*, m/s, apply 0.975 correction factor per H&H (pp. 70)
R_t = sqrt(A_t/pi); %Calculate Throat Radius (m)
R_c = R_t*sqrt(conRat); %Chamber Radius (m)
radRat = 0.7; %R2/R2max - Refer to RPA User Manual : BETWEEN 0 & 1
expRat = mainProperties(end,2); %Nozzle Expansion Ratio
conAng = 35; %Nozzle contraction angle (degrees)
plot1 = true;% Plot output in chamberContour Function call
[engineContour,chBarrel,nozzleContour, R_tCurve] = getContour(R_t,LStar,conRat,conAng,radRat,expRat,plot1);
combProperties = [mainProperties(1:2,:); nozzleProps_working];
[~,numProps] = size(combProperties);

%% Calculate Nozzle Flow Properties Across Contour
%Separate Initial Properties into Chamber and Nozzle - then later merge -
%Chamber calculations will use a linear relation from the nominal chamber to combustion
%end; Nozzle calculations will use a pchip interpolation across contraction
%ratio

%DATA FORMAT:
%pip,aeat,mach,cf,ivac,isp,p,t,rho,h,u,mw,cp,gam,son,vis,cond,condfz,pran,pranfz
AR_Nozzle = (nozzleContour(1,:).^2)./(R_t^2);
[~,throatInd] = find(nozzleContour(1,:) == R_t);
[~,throatInd1] = find(nozzleProps_working(:,2)' == 1)
[numRows,numColumns] = size(nozzleProps_working);
nozzleProps_working(1:(throatInd1-1),2) = (1-abs(1-nozzleProps_working(1:(throatInd1-1),2)));
for i = 1:numColumns %Ensure all data is unique for interpomalations
    [~,ui_new] = unique(nozzleProps_working(:,i),'stable');
    if numel(ui_new) < numRows
       numRows = numel(ui_new);
       nozzleProps_working = nozzleProps_working(ui_new,:);
    end
end
AR_Nozzle(1:(throatInd-1)) = 1-abs(1-AR_Nozzle(1:(throatInd-1)));
nozzleProps = [];
tCorr = 5; %Throat Correction Factor
for i = 1:numColumns
    newColumn = [];
    for j = 1:numel(AR_Nozzle) 
        interpUpper = interp1(nozzleProps_working(:,2)',nozzleProps_working(:,i)',AR_Nozzle(throatInd+tCorr),'linear');
        interpLower = interp1(nozzleProps_working(:,2)',nozzleProps_working(:,i)',AR_Nozzle(throatInd-tCorr),'linear');
        if abs(throatInd-j) <= tCorr %% Near throat - use special interpolation
            %So here's the deal - since M = 1 at the throat there is a
            %shock and it leads to local discontinuities. It creates flat
            %regions in the plot during interpolation as seen in
            %FigureD.1_ThroatInterpError.png. To combat this, code linearly
            %interpolates to throat set point when it's within tCorr points
            %of throat index in contour  effectively, uses Radial coordinate
            %for interpolation rather than AR
            if j < throatInd
                Z1 = [nozzleContour(2,throatInd-tCorr), nozzleContour(2,throatInd)]; %Establish interpolation AR points
                Prop1 = [interpLower, nozzleProps_working(throatInd1,i)];
                newData = interp1(Z1,Prop1,nozzleContour(2,j),'pchip');
                newColumn = [newColumn;newData];
            elseif j > throatInd
                Z1 = [nozzleContour(2,throatInd),nozzleContour(2,throatInd+tCorr)]; %Establish interpolation AR points
                Prop1 = [nozzleProps_working(throatInd1,i),interpUpper];
                newData = interp1(Z1,Prop1,nozzleContour(2,j),'pchip');
                newColumn = [newColumn;newData];
            else
                newColumn = [newColumn;nozzleProps_working(throatInd1,i)]; %Add data for throat point
            end
        else
            newData = interp1(nozzleProps_working(:,2)',nozzleProps_working(:,i)',AR_Nozzle(j),'linear');
            newColumn = [newColumn;newData];
        end
    end
    nozzleProps = [nozzleProps,newColumn];
end
nozzleProps((nozzleProps(:,2)<1),2) = abs(nozzleProps((nozzleProps(:,2)<1),2)-1)+1;
nozzleProps(end,:) = mainProperties(end,:); %Sub in exit properties for NaN results
nozzleProps = [nozzleContour(1,:)',nozzleContour(2,:)',nozzleProps]; %Add in contour vector for ideal vector
chamberProps = [chBarrel(1,:)',chBarrel(2,:)']; %Initialize vector for chamber properties
for i = 1:numColumns %Interpolate combustion properties in the chamber (#End of chamber considered to also be end of combustion - combustion process assumed linear
    newData = interp1([0,chBarrel(2,end)],[mainProperties(1,i),mainProperties(2,i)],chBarrel(2,:));%Interpolate from injector to end of chamber
    chamberProps = [chamberProps,newData'];
end
engineProps = [chamberProps;nozzleProps(2:end,:)];

%% REFERENCES %%









