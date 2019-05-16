function ndsi = NDSI( green,swir )
%NDSI calculate Normalized Difference Snow Index (NDSI) using Green and
%     SWIR bands.
% Syntax
%
%     ndsi = NDSI(green,swir)
%
% Description
%
%     This function calculates Normalized Difference Snow Index (NDSI) 
%     using Green and SWIR bands (as following equation). This can be used 
%     to detect snow, as the atmosphere is transparent at both these 
%     wavelengths, while snow is very reflective at 0.66 mm and not 
%     reflective at 1.6mm. At visible wavelengths (e.g. 0.66 microns), snow
%     cover is just as bright as clouds, and is therefore difficult to 
%     distinguish from cloud cover. However, at 1.6 microns, snow cover 
%     absorbs sunlight, and therefore appears much darker than clouds.
%
%     NDSI=(Green-SWIR)/(Green+SWIR).
%
% Input arguments
%
%     green       Green band
%     swir        Short-wave infrared band
%
% Output arguments
%
%     ndsi        Normalized Difference Snow Index
%
%
% Author: Shi Qiu (shi.qiu@ttu.edu)
% Date: 19. October, 2017

    % calculate NDSI
    ndsi=(green-swir)./(green+swir);
    
    % fix unnormal pixels
    ndsi((green+swir)==0)=0.01;

end

