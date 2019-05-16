% verison
% change to project the pixels from the potential shadow layer. Mar., 3, 2018 by Shi Qiu
% sometimes the projection will generate Nan values. fixed this bug by
% Feb., 13, 2018 by Shi Qiu
% fix the bug some pixels out of image. Shi 9/12/2017
%
function [dem_plane_r,dem_plane_c]= ProjectDEM2Plane(dim,pshadow,dem,dem_base,...
    sun_elevation_rad,sun_tazi_rad,sun_azimuth_deg,dim_expd,resolu)
    % DEM is vertically obervated by scene.
    % assume resolu.x=resolu.y
    spatial_resol=resolu(1);
    % get the location of orginal pixel  
    dem_loc=zeros(dim,'uint8')+1;
    [loc_i,loc_j]=find(dem_loc);

% must project all pixels.
%     [loc_i,loc_j]=find(pshadow==1);
%     loc_i=uint16(loc_i);
%     loc_j=uint16(loc_j);
%     if isequal(dem,0) % when no DEMs, no cloud shadow shape correction.
    if isempty(dem) % when no DEMs, no cloud shadow shape correction.
        loc_i_plane=loc_i;
        loc_j_plane=loc_j;
    else % when having DEMs
%         dem_dif=dem(pshadow==1)-dem_base;
        dem_dif=dem-dem_base;
        dem_dif=dem_dif(:);
        clear dem_loc dem dem_base;
        % sun angle geometry
%         sun_elevation_deg=90-sun_zenith_deg;
%         sun_elevation_rad=deg2rad(sun_elevation_deg);
        d_ij=dem_dif./(spatial_resol*tan(sun_elevation_rad));
%         Sun_tazi=sun_azimuth_deg-90;
%         sun_tazi_rad=deg2rad(Sun_tazi);
%         if sun_azimuth_deg< 180
        
        if sun_azimuth_deg< 180
            di=0-d_ij.*sin(sun_tazi_rad); % i for row, Y
            dj=0-d_ij.*cos(sun_tazi_rad); % j for col, X
        else
            di=d_ij.*sin(sun_tazi_rad); % i for row, Y
            dj=d_ij.*cos(sun_tazi_rad); % j for col, X
        end
        di=double(di);
        dj=double(dj);
        clear dem_dif d_ij Sun_tazi sun_elevation_rad sun_zenith_deg Sun_tazi sun_tazi_rad;

        % original location adds the bias.
        loc_i_plane=round(loc_i+di);
        loc_j_plane=round(loc_j+dj);
        clear di dj;
    end
%      dim_expd=1000;% 1000 buffer
    dim_plane_expd=dim+(2*dim_expd);
    % expand the range of observations
    loc_i_plane_expd=loc_i_plane+dim_expd;
    loc_j_plane_expd=loc_j_plane+dim_expd;

    % remove the pixels out of box.
% %     out_ids=loc_i_plane_expd(:)<1|loc_j_plane_expd(:)<1|...
% %         loc_i_plane_expd(:)>dim_plane_expd(1)|loc_j_plane_expd(:)>dim_plane_expd(2);
% %     loc_i_plane_expd(out_ids)=[];
% %     loc_j_plane_expd(out_ids)=[];
% %     loc_i(out_ids)=[];
% %     loc_j(out_ids)=[];
    
    loc_i_plane_expd(loc_i_plane_expd(:)<1)=1;
    loc_j_plane_expd(loc_j_plane_expd(:)<1)=1;
    loc_i_plane_expd(loc_i_plane_expd(:)>dim_plane_expd(1))=dim_plane_expd(1);
    loc_j_plane_expd(loc_j_plane_expd(:)>dim_plane_expd(2))=dim_plane_expd(2);
    
    clear loc_i_plane loc_j_plane dj dim_expd;
    dem_plane_r=zeros(dim_plane_expd);
    dem_plane_c=zeros(dim_plane_expd);
%     dem_plane_r = NaN(dim_plane_expd,'double');  
%     dem_plane_c = NaN(dim_plane_expd,'double'); 
    % recorders.
    tmp_id_plane=sub2ind(dim_plane_expd,double(loc_i_plane_expd),double(loc_j_plane_expd));
    clear loc_i_plane_expd loc_j_plane_expd;
    
    % sometimes there may be some Nan values!
    nan_ids=isnan(tmp_id_plane);
    if sum(nan_ids(:))>0
        tmp_id_plane(nan_ids)=[];
        loc_i(nan_ids)=[];
        loc_j(nan_ids)=[];
    end
    dem_plane_r(tmp_id_plane)=loc_i;
    dem_plane_c(tmp_id_plane)=loc_j;
    
    
    
    clear loc_i loc_j tmp_id_plane;
    % filled the holes due to the dicretization projection.
%     dem_plane_r=fillmissing(dem_plane_r,'nearest',1);% by column.
%     dem_plane_c=fillmissing(dem_plane_c,'nearest',1);% by column.
    % fill remainings as 0 value.
%     dem_plane_r=fillmissing(dem_plane_r,'constant',0);% by column.
%     dem_plane_c=fillmissing(dem_plane_c,'constant',0);% by column.

%     dem_plane_r=round(dem_plane_r);
%     dem_plane_c=round(dem_plane_c);
    
    
%     dem_plane_r=imfill(dem_plane_r,'holes');
%     dem_plane_c=imfill(dem_plane_c,'holes');
end

