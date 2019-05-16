function pfpl = DetectPotentialFalsePositivePixels(mask, psnow,slope, ndbi,ndvi,bt,cdi,water,resolution)
%DETECTPOTENTIALFALSEPOSITIVEPIXELS 
% remove the water pixels at the final layer. by Shi Mar.27, 2018.
    %%  urban
    ndbi = EnhanceLine(ndbi); % enhance the ndbi to high urban/built-up area.
    pfpl = (ndbi > 0)&(ndbi > ndvi)&mask&water==0;% urban
    
    %% remove cloudy pixels using ostu method for thermal band.
    if sum(pfpl(:))>0 % have urban pixels
        if ~isempty(bt)
            % ostu's method to determine the temperature therhold to
            % seperate clouds from urban/built-up area.
            bt_pfpl = bt(pfpl==1);
            t = graythresh(bt_pfpl);
            % get the thershold
            BW = imbinarize(bt_pfpl,t);
            clear t;
            cold_tmp = min(bt_pfpl(BW));
            clear bt_pfpl BW;
            if ~isempty(cold_tmp)
                clouds_tmp = bt < cold_tmp(1);
                pfpl(clouds_tmp==1)=0;
            end
            clear bt cold_tmp
%            figure; imshow(pfpl);
        elseif ~isempty(cdi)
%             pfpl(cdi < -1)=0; % by following David, 2018 RSE
            pfpl(cdi < -0.8)=0; % by following David, 2018 RSE
        end
    end
    
    % remove the isolate urban/built-up pixels
% %     pfpl = bwareaopen(pfpl,2);

    %% add poential snow/ice pixels in mountain.
    if ~isempty(slope)
        psnow_mountain = psnow==1& slope > 20;
        % 20 is from Burbank D W, Leland J, Fielding E, et al. Bedrock incision, rock uplift and threshold hillslopes in the northwestern Himalayas[J]. Nature, 1996, 379(6565): 505.
        pfpl = pfpl|psnow_mountain;
    end
    
    %% buffer urban pixels with 1 kilometer window. 
    witdh=250;
    width_pixels=fix(witdh/resolution);% 1 km 33 Landsat pixel   500 meters 17 Landsat pixels.   200 meters 7 pixels.
    SEw=strel('square',2*width_pixels+1);
    pfpl=imdilate(pfpl,SEw);
    
    %% must over land this will be runned in commission process function.
    % some coastline may be classified as water, that is easily detected as
    % cloud when it is grouped into water layer.
%     pfpl = pfpl&mask&water==0; % remove the water pixels at the final. also remove the absolute cloud pixels.
    
% %     % dilate user-defined buffers for erosion more commission errors.
% %     if pfbpix > 0
% %        SEw=strel('square',2*pfbpix+1);
% %        pfpl=imdilate(pfpl,SEw);
% %     end
    
    pfpl = pfpl|psnow;%% add the snow/ice pixels located in normal areas.
    
    pfpl = pfpl&mask; % remove the water pixels at the final. also remove the absolute cloud pixels.
end