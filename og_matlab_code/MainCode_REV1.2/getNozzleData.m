function [mainProperties,nozzleProps] = getNozzleData(P_rat2,conRat,CEAData_int,pRatStations,areaRatStations)
%getNozzleData: Generates combustion property data across nozzle based on
%Area Ratio
%   Detailed explanation goes here
% Outputs
%   - Main Properties (Chamber, Comb End, Throat,Nozzle Exit)
%   - Converging Properties (Comb End to Throat based on Ae/At)
%   - Diverging Properties (Throat to Nozzle Exit based on Ae/At)
%
% DATA LAYOUT: 
%pip,aeat,mach,cf,ivac,isp,p,t,rho,h,u,mw,cp,gam,son,vis,cond,condfz,pran,pranfz
CEAData_int = CEAData_int(:,2:end); %Strip constant CR Measurement
mainProperties = CEAData_int(1:3,:); %Does not yet include Nozzle Exit
convergingProperties_raw = [mainProperties(2,:);CEAData_int(4:(3+numel(pRatStations)),:)];
divergingProperties_raw = [CEAData_int((end+1-numel(areaRatStations)):end,:)];

% CONVERGING PROPERTIES:     
conv_Arat = linspace(conRat,1,100); %Query Area Ratios for Interpolation
convergingProperties = []; %Final Converging Props Array
minRow = find(mainProperties(3,1) == sort([mainProperties(3,1),convergingProperties_raw(:,1)']))-1;
convergingProperties_raw = [convergingProperties_raw(1:minRow,:);mainProperties(3,:)]; %Format data to end at throat and cut off data before throat
interp_Arat1 = convergingProperties_raw(:,2); %'X' Coordinate for Interpolation
[numRows,numColumns] = size(convergingProperties_raw);
for i = 1:numColumns
    newData = interp1(interp_Arat1',convergingProperties_raw(:,i)',conv_Arat,'pchip');
    convergingProperties = [convergingProperties,newData'];
end

% DIVERGING PROPERTIES
div_Arat = linspace(1,divergingProperties_raw(end,2),100); %Query Area Ratios for Interpolation
divergingProperties = []; %Final Diverging Props Array
[numRows,numColumns] = size(convergingProperties_raw);
divergingProperties_raw = [mainProperties(3,:);divergingProperties_raw(:,:)];
interp_Arat2 = divergingProperties_raw(:,2); %'X' Coordinate for Interpolation
for i = 1:numColumns
    newData = interp1(interp_Arat2',divergingProperties_raw(:,i)',div_Arat,'pchip');
    divergingProperties = [divergingProperties,newData'];
end
endRow = find(P_rat2 == sort([P_rat2,divergingProperties(:,1)']))-1; %Find index of last data before exit
expRat = interp1(divergingProperties(endRow:(endRow+1),1)',divergingProperties(endRow:(endRow+1),2)',P_rat2); %Solve for expansion ratio based on pressure
exitProperties = []; %Vector of exit conditions for combustion gasses
for i = 1:numColumns %Interpolate via Area ratio for other properties to maintain consistency.
    newData = interp1(divergingProperties(endRow:endRow+1,2)',divergingProperties(endRow:endRow+1,i)',expRat);
    exitProperties(i) = newData;
end
divergingProperties = [divergingProperties(1:endRow,:);exitProperties]; %Add exit properties as final row
mainProperties = [mainProperties;exitProperties];
mainProperties(1,2) = conRat;

%% R&D
convergingProperties_rd = convergingProperties;
convergingProperties_rd(:,2) = (1-abs(1-convergingProperties_rd(:,2)));
nozzleProps_raw = [convergingProperties_rd(1:(end-1),:);divergingProperties]; %Combine converging and diverging properties
uniqueIndicies=[];
[numRows,numColumns] = size(nozzleProps_raw);
for i = 1:numColumns %Ensure all data is unique for interpomalations
    [~,ui_new] = unique(nozzleProps_raw(:,i),'stable');
    if numel(ui_new) < numRows
       numRows = numel(ui_new);
       nozzleProps_raw = nozzleProps_raw(ui_new,:);
    end
end
conAR = (1-abs(1-linspace(conRat,1,100)));
divAR = linspace(1,expRat,101);
ARQueries = [conAR,divAR(2:end)];
nozzleProps = [];
for i = 1:numColumns
    newData = interp1(nozzleProps_raw(:,2)',nozzleProps_raw(:,i)',ARQueries,'pchip');
    %newData = ppval(ppol,ARQueries);
    nozzleProps = [nozzleProps,newData'];
end
nozzleProps((nozzleProps(:,2)<1),2) = abs(nozzleProps((nozzleProps(:,2)<1),2)-1)+1;
% plot(nozzleProps(:,2),nozzleProps(:,8),'.r')
% 'stop';

end

