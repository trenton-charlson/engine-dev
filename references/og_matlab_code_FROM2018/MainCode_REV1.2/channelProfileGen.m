function [channelProfile] = channelProfileGen(inputFile)
%channelProfileGen - generates matrix of relevant channel geometric
%properties from input excel spreasheet

%Design inputs: 
% - Chamber Wall Thickness
% - Fin width
% - Num Channels
% - R, Z

%Design Outputs
% - Channel CSA
% - Channel Arc Len

%Spreadsheet Layout:
% R, Z, chamberWall, channelFloor, floorCircum, channelWidth, numChannels,
% finWidth, finHeight (depth), channelCSA, channelAR, channelPerim, outerJacket
%Note - CSA is currently rough estimate
channelProfile_inp = xlsread(inputFile);

chamberRad = channelProfile_inp(:,1);
channelDepth = channelProfile_inp(:,9);
channelFloorRad = channelProfile_inp(:,4);
channelTopRad = channelProfile_inp(:,13);
channelWidth = channelProfile_inp(:,6);
numChannels = channelProfile_inp(:,7);
finWidth = channelProfile_inp(:,8);
channelAR = channelProfile_inp(:,11);
channelCSA = channelProfile_inp(:,10);
channelPerim = channelProfile_inp(:,12);

%Output Layout
% R, Z, R_floor, R_ceil, numChannels, channelWidth, finWidth, channelDepth, 
%channelCSA, channelPerim, channelAR
channelProfile = [chamberRad, channelProfile_inp(:,2), channelFloorRad, channelTopRad, numChannels, channelWidth, finWidth, channelDepth, channelCSA, channelPerim, channelAR];

end

