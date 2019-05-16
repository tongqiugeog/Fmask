function ndvi = NDVI( red,nir )
%NDVI Calculate Normalized Difference Vegetation Index (NDVI) using NIR and
%     Red bands.
%
% Syntax
%
%     ndvi = NDVI(red,nir)
%
% Description
%
%     This function calculates Normalized Difference Vegetation Index (NDVI) 
%     using NIR and Red bands (as following equation). This range is between
%     -1 and 1.
%     NDVI=(NIR-Red)/(NIR+Red).
%
% Input arguments
%
%     red         Red band
%     nir         Near-infrared band
%
% Output arguments
%
%     ndvi        Normalized Difference Vegetation Index
%
%
% Author: Shi Qiu (shi.qiu@ttu.edu)
% Date: 19. October, 2017

    % calculate NDVI
    ndvi=(nir-red)./(nir+red);
    
    % fix unnormal pixels
    ndvi((nir+red)==0)=0.01;

end

