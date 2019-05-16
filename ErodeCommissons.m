
% output: cloud 1:cloud; 2: potential cloud can be eroded, which will be
% identified using cloud shadow match method.
function cloud = ErodeCommissons(data_meta,pcloud,pfpl,water,cdi,erdpix)
%REMOVECOMMISSONS remove most of commission errors from bright and low
%tempareture features.
% 
% use the optimal erosion size for Landsat and Sentinel 2 images. Shi
% 4/21/2018.
% remove the cloud objects with less than 3 pixels after eroding 1 pixel.
% Shi 4/10/2018
 % pcloud potential cloud layer
 % pccl potential commissions as cloud layer

    % 2-by-2 pixels window  for Landsat 8
    % 3-by-3 pixels window  for Landsats 4-7 and Sentinel 2
    cipix = erdpix;
    
    cloud = zeros(data_meta.Dim,'uint8');  % the clouds needed to be eroded.
    cloud(pcloud>0)=1;
    clear pcloud;
    
    
    %% erode and dilate back.
%     CEs = strel('square',2*cipix+1);
    CEs = strel('disk',cipix);
    % erode to remove the potential false pixels
    cloud_eroded = imerode(cloud,CEs);
    clear CEs;
    pixels_eroded = ~cloud_eroded & pfpl;
    clear cloud_eroded;
    % only remove the potential false positive pixels
    cloud_eroded = cloud;
    cloud_eroded(pixels_eroded) = 0;
    clear pixels_eroded;
    
    %% dilate back to orginal cloud shape of which size is dual to the erosion.
    CEs = strel('disk',2*cipix);
    cloud_dilated = imdilate(cloud_eroded,CEs);
    
    %% remover the clouds gone forever.
    % Segmentate each cloud to remove the small objs.
    cloud_objs=bwlabeln(cloud,8);
    clouds_remaining = cloud_objs;
    clouds_remaining(~cloud_eroded)=0; % remove the clouds gone.
    clear cloud_eroded;
    idx_clouds_remaining = unique(clouds_remaining(:));
    clear clouds_remaining;
    
    cloud_remaining = zeros(data_meta.Dim,'uint8');  % the clouds needed to be eroded.
    if ~isempty(idx_clouds_remaining)
        if idx_clouds_remaining(1)==0
            idx_clouds_remaining(1) = [];
        end
        if ~isempty(idx_clouds_remaining)
            cloud_remaining = ismember(cloud_objs,idx_clouds_remaining);
        end
    end
    clear cloud_objs idx_clouds_remaining;
    
    % only for land
    
    clear cloud_eroded;
    cloud = (cloud_dilated&cloud_remaining)|(water&cloud); % add clouds over water.
	clear cloud_dilate cipix CEs water;
    
    %% remove small object with minum CDI < 0.5 only for Sentinel 2
    if ~isempty(cdi)
        % exclude small objects of which minimum CDI is still larger than
        % -0.5.
        % get small clouds
        large_obj = bwareaopen(cloud,10000);
        small_obj = cloud==1&large_obj==0;
        clear large_obj;
        % segment small clouds
        small_obj_init=bwlabeln(small_obj,8);
        
        % produce all non confident cloud pixels using CDI
        confident_cloud = cdi < -0.5;
        
        clear cdi;
        
        % remove the non confident cloud pixels again
        small_obj_exd_urban = small_obj_init;
        small_obj_exd_urban(confident_cloud==0) = 0;
        clear confident_cloud;
        
        % a true cloud can be determined when any confident pixels are
        % remaining.
        idx = unique(small_obj_exd_urban);
        clear small_obj_exd_urban;
        
        true_cloud = ismember(small_obj_init,idx);
        clear small_obj_init idx;
        
        % remove bright surfaces
        bright_surfaces = true_cloud==0&small_obj==1;
        clear small_obj true_cloud;
        cloud(bright_surfaces) =0;
    end
    
    %% remove very small object.
    cloud = bwareaopen(cloud,3);
    
end