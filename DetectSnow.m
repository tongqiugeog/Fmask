function snow = DetectSnow(dim,band_green,band_nir,band_bt,ndsi)
%DETECSNOW Summary of this function goes here
%   Detailed explanation goes here
  
% %     snow=zeros(dim,'uint8'); % Snow mask
    % It takes every snow pixels including snow pixel under thin clouds or icy clouds
    snow=ndsi>0.15&band_nir>1100&band_green>1000;
    if ~isempty(band_bt)
        snow=snow&band_bt<1000;
    end
%     snow(ids_snow)=1;
end

