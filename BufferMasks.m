function [cloud,shadow,snow] = BufferMasks(pcloud,cdpix,pshadow,csdpix,psnow,sdpix)
% BUFFERMASKS Dilate cloud/cloud shadow/snow with cdpix/csdpix/sdpix with
%     default value of 3/3/3.
%
% Syntax
%
%     [cloud,shadow,snow] =
%     BufferMasks(pcloud,cdpix,pshadow,csdpix,psnow,sdpix)
%
% Description
%
%     Dilate cloud/cloud shadow/snow with cdpix/csdpix/sdpix with
%     default value of 3/3/3.
%
% Input arguments
%     pcloud    Cloud mask.
%     cldpix    Dilated number of pixels for cloud with default value of 3.
%     pshadow   Cloud shadow mask.
%     sdpix     Dilated number of pixels for cloud shadow value of 3.
%     psnow     Snow mask.
%     snpix     Dilated number of pixels for snow value of 3.
%
% Output arguments
%
%     cloud    Final dilated cloud mask.
%     shadow   Final dilated cloud shadow mask.
%     snow     Final dilated snow mask.
%
% Example
%
%     [cloud,shadow,snow] =
%     BufferMasks(pcloud,3,pshadow,3,psnow,3);
%
%        
% Author:  Shi Qiu (shi.qiu@ttu.edu)
% Date: 2. November, 2017 

    % buffer cloud
    if cdpix>0
        CEs=strel('square',2*cdpix+1);
        cloud=imdilate(pcloud,CEs);
        clear pcloud CEs cdpix;
    else
        cloud=pcloud;
        clear pcloud
    end
    % buffer cloud shadow
    if csdpix>0
        CSEs=strel('square',2*csdpix+1);
        shadow=imdilate(pshadow,CSEs);
        clear pshadow CSEs csdpix;
    else
        shadow=pshadow;
        clear pshadow
    end
    % buffer snow
    if sdpix>0
        SEs=strel('square',2*sdpix+1);
        snow=imdilate(psnow,SEs);
        clear psnow SEs sdpix;
    else
        snow=psnow;
        clear psnow
    end
end

