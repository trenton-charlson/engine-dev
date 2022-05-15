function [CEAData_int] = getDataDeck(dataMain,P_c,conRat,conRat_array,pStations)
%getDataDeck: Import, retrieve, and interpolate data from given CEA Data
%Deck specified by user

dataIndex = find(conRat == conRat_array);
if numel(dataIndex) > 0 %Contraction Ratio is one queried
    dataQuery = dataMain{dataIndex};
    [totalLen, ~] = size(dataQuery);
    unitLength = totalLen/numel(pStations);
    if any(P_c == pStations) %Desired P_c in data deck queries
        %Grab data for given Pc (and CR - baseline case)
        pCi = find(P_c == pStations);
        CEAData_int = dataQuery(((pCi-1)*unitLength)+1:(pCi*unitLength),1:end);
    else
        %Grab low and high data for interpolation
        pCi_high = find(P_c == sort([P_c pStations]));
        pCi_low = pCi_high-1;
        dataLow = dataQuery(((pCi_low-1)*unitLength)+1:(pCi_low*unitLength),:);
        dataHigh = dataQuery(((pCi_high-1)*unitLength)+1:(pCi_high*unitLength),:);
        CEAData_int = [];
        [~,numColumns] = size(dataLow);
        pC_Data = [dataLow(1,8),dataHigh(1,8)];
        for i = 1:unitLength %Interpolate between chamber pressures
            for j = 1:numColumns
                dataQuery = [dataLow(i,j),dataHigh(i,j)];
                CEAData_int(i,j) = interp1(pC_Data,dataQuery,P_c,'pchip');
            end
        end
    end        
else %Interpolation in between contraction ratio must be done
    cRi_high = find(conRat == sort([conRat, conRat_array]));
    cRi_low = cRi_high-1;
    data_cRi_low = dataMain{cRi_low};
    data_cRi_high = dataMain{cRi_high};
    cR_Data = [conRat_array(cRi_low), conRat_array(cRi_high)];
    [totalLen, ~] = size(data_cRi_low);
    unitLength = totalLen/numel(pStations);
    if any(P_c == pStations) %Desired P_c in data deck queries
        %Grab data for given Pc from each cR
        pCi = find(P_c == pStations);
        dataLow = data_cRi_low(((pCi-1)*unitLength)+1:(pCi*unitLength),1:end);
        dataHigh = data_cRi_high(((pCi-1)*unitLength)+1:(pCi*unitLength),1:end);
        CEAData_int = [];
        [~,numColumns] = size(dataLow);
        for i = 1:unitLength % Interpolate between contraction Ratios
            for j = 1:numColumns
                dataQuery = [dataLow(i,j),dataHigh(i,j)];
                CEAData_int(i,j) = interp1(cR_Data,dataQuery,conRat,'pchip');
            end
        end 
    else
        %Interpolate between Pc and cR
        pCi_high = find(P_c == sort([P_c pStations])); %index of higher chamber pressure in data deck
        pCi_low = pCi_high-1;
        dataLow1 = data_cRi_low(((pCi_low-1)*unitLength)+1:(pCi_low*unitLength),:); %For lower cR
        dataHigh1 = data_cRi_low(((pCi_high-1)*unitLength)+1:(pCi_high*unitLength),:);
        dataLow2 = data_cRi_high(((pCi_low-1)*unitLength)+1:(pCi_low*unitLength),:); %for upper cR
        dataHigh2 = data_cRi_high(((pCi_high-1)*unitLength)+1:(pCi_high*unitLength),:);
        [~,numColumns] = size(dataLow1);
        data_int1 = []; %Intermediate data 1
        data_int2 = [];
        CEAData_int = [];
        pC_Data = [dataLow1(1,8),dataHigh1(1,8)];
        for i = 1:unitLength %Interpolate at lower contraction ratio
            for j = 1:numColumns
                dataQuery = [dataLow1(i,j),dataHigh1(i,j)];
                data_int1(i,j) = interp1(pC_Data,dataQuery,P_c,'pchip');
            end
        end
        for i = 1:unitLength %Interpolate at upper contraction ratio
            for j = 1:numColumns
                dataQuery = [dataLow2(i,j),dataHigh2(i,j)];
                data_int2(i,j) = interp1(pC_Data,dataQuery,P_c,'pchip');
            end
        end
        for i = 1:unitLength %Interpolate between lower and upper contraction ratio with intermediate data
            for j = 1:numColumns
                dataQuery = [data_int1(i,j),data_int2(i,j)];
                CEAData_int(i,j) = interp1(cR_Data,dataQuery,conRat,'pchip');
            end
        end
        
    end 
end

end

