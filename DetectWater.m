function [water, waterAll] = DetectWater( dim, mask, nir, NDVI, psnow, slope, gswater)
%DETECTWATER Detect water by combining spectral-derived water and
%GSWO-derived water togeter.
%
% Syntax
%
%     water = DetectWater( dim, mask, nir, NDVI, psnow, slope, gswater)
%
% Description
%
%     History:
%     1. Create this function. (1. January, 2018)
%     2. The sepctral-derived water may be incorrect, resulting in a 100%
%     absolutely wrong low level GSWO (equal to 0). The GWSO will be used
%     only when the low level GSWO is larger than 0. (9. March, 2018)
%     3. Remove the coastline because of its frequent changes.  (6. May, 2018)
%     4. Add a water layer which does not snoe/ice because some clouds may
%     be like snoe/ice. This will be used to exclude processing of cloud
%     shadow over water. (17. March, 2020)
%
%
% Input arguments
%
%     dim            Dim for data.
%     mask           Mask for observations.
%     nir            NIR.
%     NDVI           NDVI.
%     psnow          Potential snow.
%     slope          Slope.
%     gswater        GSWO.
%
% Output arguments
%
%     water          Water map.
%
%        
% Author:  Shi Qiu (shi.qiu@uconn.edu)
% Date: 17. March, 2020
  
    water=zeros(dim,'uint8'); % Water msk
    %% Zhe's water test (works over thin cloud)
    water((NDVI<0.01&nir<1100)|(NDVI<0.1&NDVI>0&nir<500))=1;
    clear resolution;
    
    % within observation.
    water(~mask)=0; 
    % do not exclude snow over water because clouds may be like snow and
    % will be excluded ...
    waterAll = water;
    %% the GSWO data to enhance water map.
    if sum(water(:))>0&&~isempty(gswater)
        if sum(gswater(:))>0 % have water occurences.
            % assume the water occurances are same in each whole scene.
            % global surface water occurance (GSWO)
            % low level to exclude the commssion errors as water.
            % 5% tolerances.
            gswater_occur=prctile(gswater(water==1),17.5)-5;
            
            if gswater_occur>0 % must be more than 0.
                water_gs = gswater>gswater_occur;
                clear gswater gswater_occur;
                
                waterAll(water_gs==1)=1;% extend water regions based on GSWO, but do not exclude snow/ice
                
                % sometimes ice may be over water. Snow covered sea ice is determined
                % using the potential snow/ice.
                water_gs(psnow == 1)=0; % remove out ice water.
                clear psnow;
                
                water(water_gs==1)=1;% extend water regions based on GSWO. 
                % equal to the theshold because it may be 100%.
                
    %             water(psnow)=0; % remove out ice water. I think this snow
    %             cannot be removed because there are sometimes ice clouds over
    %             water.
                water(~mask)=0;
                waterAll(~mask)=0;
            end
        end
        % note that 255 indicates no data in GSWO, that is ocean pixels or
        % permenate snow/ice pixels (which can be identified as land pixels).
    end

end

