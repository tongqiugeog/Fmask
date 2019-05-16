function ndbi = NDBI( nir, swir )
%NDBI Calculate Normalized Difference Build-up Index (NDBI) using NIR and
%     SWIR bands.
% Syntax
%
%     ndbi = NDBI(swir,nir)
%
% Description
%
%     This function calculates Normalized Difference Build-up Index (NDBI) 
%     using SWIR and NIR bands (as following equation).
%
% Input arguments
%
%     swir        Short-wave infrared band
%     nir         Near infrared band
%
% Output arguments
%
%     ndbi        Normalized Difference Build-up Index
%
%
% Author: Shi Qiu (shi.qiu@ttu.edu)
% Date: 19. Dec., 2017

    % calculate NDBI
    ndbi=(swir-nir)./(swir+nir);
    
    % fix unnormal pixels
    % not 0.01 any more because we will identify urban pixel using ndbi more than 0.
%     ndbi((nir+swir)==0)=0.0;

end

