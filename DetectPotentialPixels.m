function [idplcd,BandCirrusNormal,whiteness,HOT] = DetectPotentialPixels(mask,data_toabt,dem,ndvi,ndsi,satu_Bv)
% DETECTPOTENTIALCLOUD detect potential cloud pixels (PCPs)

    % Cirrus Probability  This is unavailable here because some high
    % mountianus have high cirrus values.
    
    % inputs: BandSWIR2 BandBT BandBlue BandGreen BandRed BandNIR BandSWIR1
    % BandCirrus
    %% Step 1: detect possible cloud pixels (PCPs)
    
    % [Basic cloud test]
    idplcd=ndsi<0.8&ndvi<0.8&data_toabt.BandSWIR2>300;    
    clear ndsi ndvi; data_toabt.BandSWIR2 = []; % memory.
    % when have BT data.
    if ~isempty(data_toabt.BandBT)
        idplcd=idplcd==true&data_toabt.BandBT<2700;
        data_toabt.BandBT = [];
    end
    
    % [Whiteness test]
    % visible bands flatness (sum(abs)/mean < 0.6 => brigt and dark cloud )
    visimean=(data_toabt.BandBlue+data_toabt.BandGreen+data_toabt.BandRed)/3;
    whiteness=(abs(data_toabt.BandBlue-visimean)+abs(data_toabt.BandGreen-visimean)+...
        abs(data_toabt.BandRed-visimean))./visimean;
    data_toabt.BandGreen = [];
    clear visimean;
    % update idplcd
    whiteness(satu_Bv==1)=0;% If one visible is saturated whiteness == 0
    idplcd=idplcd==true&whiteness<0.7;

    % [Haze test]
    HOT=data_toabt.BandBlue-0.5.*data_toabt.BandRed-800;% Haze test
    data_toabt.BandBlue = [];
    data_toabt.BandRed = [];

    idplcd=idplcd==true&(HOT>0|satu_Bv==1);
    clear satu_Bv; % need to find thick warm cloud
    
    % [Ratio4/5>0.75 cloud test]
    Ratio4_5=data_toabt.BandNIR./data_toabt.BandSWIR1;
    data_toabt.BandNIR = [];
    data_toabt.BandSWIR1 = [];
    idplcd=idplcd==true&Ratio4_5>0.75;
    clear Ratio4_5;
    
    BandCirrusNormal=[];
    % normalize Cirrus band [Cirrus test] from Landsat 8 and Sentinel 2 images
    if ~isempty(data_toabt.BandCirrus)
        BandCirrusNormal=NormalizaCirrusDEM( mask, idplcd, data_toabt.BandCirrus, dem );
%         BandCirrusNormal= data_toabt.BandCirrus;
        clear data_toabt mask dem;
        idplcd=idplcd==true|BandCirrusNormal > 100; % When TOA at Cirrus band is more than 0.01, it may be cloudy.
    end
    
end

