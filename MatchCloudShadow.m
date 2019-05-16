%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This functin is used to match cloud shadow with cloud.
% The similarity defineing cloud shadow is larger than 0.3 only for those
% clouds, of which all pixels are included in the Landsat observations.
% This minor modification was made because the match similarity may be
% wrong when some parts of cloud are out of the observations.
% 
%
% fix the bug that cloud shadow would be projected on the other side in Sentinel-2 imagery when the azimuth angle > 180. By Shi, at 19, Jan., 2019
% use new match similarity becasue we do not know the potential clouds 
% excluding self cloud and outsides are shadow or not.   by Shi, at 21, April, 2018
% remove the overlap between final matched cloud shadow and the potential
% cloud shadow. by Shi, at 26, Mar., 2018.
% speed up the match of cloud shadow with cloud for large clouds using
% sampling projections. by Zhe and Shi. at 24, Mar., 2018.
% do not revisit for the big clouds (more than 10,000,000). by Shi. at 22, Mar., 2018
% match cloud shadow by following the sort from the center because the
% clouds loacted boundary will be easily affected due to the unkown of the
% pixels out of obervations.   by Shi. at 15, Mar., 2018
% cloud's temperature may be warmer than surface when we wrongly give some
% surface to the cloud. This will result in no cloud shadow.    by Shi. at 3, Mar., 2018
% fix the bugger that revisit clouds when less than 14 cloud objects.   by Shi. at 11, Dec., 2017
% still improve the prediction of cloud shadow location when no DEMs.  by Shi. at 13, Sept., 2017
% revisit the first 14 cloud objects.   by Shi. at 21,Feb.,2017
% fix the bugger, struct2table for lt. struct2table. at 21,Feb.,2017
% search neighboring clouds by distance rule. by Shi. at 13,Feb.,2017
% fix the bugger that bias for the location of real cloud object. at 13,Oct.,2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ similar_num,data_cloud_matched, data_shadow_matched] = MatchCloudShadow(...
    mask,plcim,plsim,pfpl,water,data_dem,data_bt_c,t_templ,t_temph,data_meta,ptm,num_near,angles_view)
    
    dim=data_meta.Dim;
    % get potential mask values
	data_shadow_potential=zeros(dim,'uint8');
    data_cloud_potential=(plcim>0)&mask==1;
    data_shadow_potential(plsim==1)=1;% plshadow layer
    clear plsim; % empty memory
    % matched cloud & shadow layer
    data_cloud_matched=zeros(dim,'uint8');
    data_shadow_matched=zeros(dim,'double');
    % revised percent of cloud on the scene after plcloud
    revised_ptm=sum(data_cloud_potential(:))/sum(mask(:));
    % When having too many clouds, we will not match cloud shadow.
    if ptm <=40000||revised_ptm>=0.90 % 0.1% was changed to 40,000 pixels.
        fprintf('Skip cloud shadow detection because high cloud cover\n');
        data_cloud_matched(data_cloud_potential==true)=1;
        data_shadow_matched(data_shadow_potential==false)=1;
        data_shadow_matched=uint8(data_shadow_matched);
        similar_num=-1;
    else
        
        clear pfpl;
        %% parameters
        clear plcim; % empty memory
        % max similarity threshold
        max_similar = 0.95;
        % number of inward pixes (240m) for cloud base temperature
        num_pix=3; 
        
        % enviromental lapse rate 6.5 degrees/km
        rate_elapse=6.5;
        % dry adiabatic lapse rate 9.8 degrees/km
        rate_dlapse=9.8;
        
        % max match pixels number
        max_match_num =1000000; %more than 1 million pixels will result in ~2 mins runtime.
        
        % sun angles
        sun_zenith_deg=data_meta.Zen;
        sun_azimuth_deg=data_meta.Azi;
        % sun angle geometry
        sun_elevation_deg=90-sun_zenith_deg;
        sun_elevation_rad=deg2rad(sun_elevation_deg);
        % solar azimuth anngle
        Sun_tazi=sun_azimuth_deg-90;
        sun_tazi_rad=deg2rad(Sun_tazi);
        clear sun_elevation_deg sun_elevation_deg;

        % view angles for Sentinel 2 images, which will be used compute the 
        % average values for each cloud. Note that the Landsat's view angles 
        % can be estimated by the obersations of the entire scene.
        if strcmp(data_meta.Sensor,'S_MSI')
            VAA = angles_view.VAA;
            VZA = angles_view.VZA ;
            clear angles_view;
            % mini matched similarity
            Tsimilar=0.425;
%             Tsimilar=0.4;
            % threshold for matching buffering
            Tbuffer=0.90; 
        elseif strcmp(data_meta.Sensor,'L_OLI_TIRS')||...
                strcmp(data_meta.Sensor,'L_ETM_PLUS')||...
                strcmp(data_meta.Sensor,'L_TM')
            % mini matched similarity
%             Tsimilar=0.35;
            Tsimilar=0.3;
            % threshold for matching buffering
            Tbuffer=0.95; 
            
            % view angle geometry for Landsat
             % get track derection
            [rows,cols]=find(mask==1);
            [cor_ul_y,num]=min(rows);cor_ul_x=cols(num);
            [cor_lr_y,num]=max(rows);cor_lr_x=cols(num);
            [cor_ll_x,num]=min(cols);cor_ll_y=rows(num);
            [cor_ur_x,num]=max(cols);cor_ur_y=rows(num);
            % get view angle geometry
            [A,B,C,omiga_par,omiga_per]=getSensorViewGeo(cor_ul_x,cor_ul_y,cor_ur_x,cor_ur_y,cor_ll_x,cor_ll_y,cor_lr_x,cor_lr_y);
            clear cor_ul_x cor_ul_y cor_ur_x cor_ur_y cor_ll_x cor_ll_y cor_lr_x cor_lr_y;
        else
            error('Only Landsats 4-7, Landsat 8 and Sentinel 2 data can be supported./n');
        end
        
        % the lowest elevation.
        if ~isempty(data_dem)
            dem_base_heigh=double(prctile(data_dem(mask),0.001));
        else
            dem_base_heigh=0;
        end
        % expand 1,000 pixels for the potential cloud shadow layer.
        dim_expd=2000; 

        % spatial resolution of the image
        sub_size=data_meta.Resolution(1);
        win_height=dim(1);win_width=dim(2);
        % intervals within each matching process.
        step_interval=2*sub_size*tan(sun_elevation_rad);

        %% project all potential cloud shadow and cloud (can be matched) along sun light based on DEMs.
        % cloud shadow may be overlap with another cloud, so we need to 
        % project the all potential cloud and potential cloud shadow pixels.
        [recorderRow,recorderCol] = ProjectDEM2Plane(dim,...
            mask,...
            data_dem,dem_base_heigh,sun_elevation_rad,sun_tazi_rad,...
            sun_azimuth_deg,dim_expd,...
            data_meta.Resolution);
        
        % create cloud objtects using 8-by-8 pixels connection.
        [segm_cloud,num]=bwlabeln(data_cloud_potential,8);
        s = regionprops(segm_cloud,'area');
        area_final = [s.Area];
        obj_num=area_final;
        clear segm_cloud_init L idx area_final s;
        
        % Get the x,y of each cloud
        % Matrix used in recording the x,y
        stats= regionprops(segm_cloud,'Centroid','PixelList');
  
        
        match_clds=zeros(1,num,'uint8'); % cloud shadow match similarity (m)
        matched_clds_centroid=[]; % centers of cloud having shadow
        height_clds_recorder=[]; % cloud shadow match heights (m)
        % Use iteration to get the optimal move distance
        % Calulate the moving cloud shadow 
        similar_num=zeros(1,num); % cloud shadow match similarity (m)
        l_pt=0.175; h_pt=1-l_pt;
        dim_expand=dim+2*dim_expd;
        record_base_h_num=0;
        num_revisited = 0;
        if num > num_near
            num_revisited=num_near;
        end
        num_all=num+num_revisited;
        
        % min moving distance (min high 200 meters) unit: pixels
        i_xy_min=200/(sub_size*tan(sun_elevation_rad)); 
        
        for cloud_type_cur= 1:num_all %num
            % revisit the first 14 cloud objects.
            cloud_type=cloud_type_cur;
            if cloud_type>num && num_revisited<num
                cloud_type=cloud_type_cur - num; % fix this bug.
            end
            
            % moving cloud xys
            XY_type_all=zeros(obj_num(cloud_type),2);
%             % record the max threshold moving cloud xys
%             tmp_XY_type_all=zeros(obj_num(cloud_type),2);
            % corrected for view angle xys
            tmp_xys_all=zeros(obj_num(cloud_type),2);
            % record the original xys
            orin_xys_all=zeros(obj_num(cloud_type),2);
            % record the original xys
            orin_xys_all(:,:)=stats(cloud_type,:).PixelList(:,:);
            % record this orinal ids
            orin_cid_all=sub2ind(dim,orin_xys_all(:,2),orin_xys_all(:,1)); 

            % assume object is round r_obj is radium of object 
            r_obj=sqrt(obj_num(cloud_type)/2*pi);

            % refine cloud height range (m)
            % initialize height and similarity info
            if isempty(data_dem)
                base_heigh_cloud=0;
            else
                % if the above rule removed all pixels, back to MFmask.
                dem_base_cloud=data_dem(orin_cid_all);    
                % The min height should be the max dem of dem_base_cloud.
                % However, the commission error from snow may lead to overestimate the base heigh. 
                base_heigh_cloud=prctile(dem_base_cloud,100*h_pt)-dem_base_heigh;
                clear dem_base_cloud;
            end
            % Min cloud base height (m)
            Min_cl_height=200.00 + base_heigh_cloud; % 2738
            % Max cloud base height (m)
            Max_cl_height=12000.00 + base_heigh_cloud;

            if ~isempty(data_bt_c) % if have no temperature data.
                % Temperature of the cloud object
                temp_obj_all=data_bt_c(orin_cid_all);
                % the base temperature for cloud
                % number of inward pixes for correct temperature
            %        num_pix=8;
                pct_obj=(r_obj-num_pix)^2/r_obj^2;
                pct_obj=min(pct_obj,1); % pct of edge pixel should be less than 1
                t_obj=quantile(temp_obj_all(:),pct_obj); 
                clear pct_obj;
                t_obj=double(t_obj);
                % put the edge of the cloud the same value as t_obj
                temp_obj_all(temp_obj_all>t_obj)=t_obj;
                if ~(t_templ<t_obj||t_temph<t_obj) % cloud's temperature may be warmer than surface when we wrongly give some surface to the cloud.
                    Min_cl_height=max(Min_cl_height,10*(t_templ-400-t_obj)/rate_dlapse);
                    Max_cl_height=min(Max_cl_height,10*(t_temph+400-t_obj));
                end
            end

            % when reaching big clouds, the max bias for cloud shadow will
            % be estimated, but exclude dem info.
            
            if obj_num(cloud_type) > max_match_num 
                % renew the arrays
                % randomly selection.
                samples_rand_all=randperm(obj_num(cloud_type)); 
                samples_mov=samples_rand_all(1:max_match_num);
                clear samples_rand_all;
                % moving cloud xys
                XY_type=XY_type_all(samples_mov,:);
                % corrected for view angle xys
                tmp_xys=tmp_xys_all(samples_mov,:);
                % record this orinal ids
                orin_xys = orin_xys_all(samples_mov,:);
%                 orin_xys = orin_xys_all(samples_mov,1);
                if ~isempty(data_bt_c) % if have no temperature data.
                    % Temperature of the cloud object
                    temp_obj=temp_obj_all(samples_mov);
                end
                
            else
                % give all pixels
                % moving cloud xys
                XY_type=XY_type_all;
                % corrected for view angle xys
                tmp_xys=tmp_xys_all;
                % record the original xys
                orin_xys=orin_xys_all;
                % record this orinal ids
%                 orin_cid=orin_cid_all; 
                
                if ~isempty(data_bt_c) % if have no temperature data.
                    % temperature                
                    temp_obj=temp_obj_all;
                end
            end

%             record_h=0;
            record_thresh=0;
            record_base_h=0;
            record_base_h_near=0;% it is available only when >0 
            center_cur=stats(cloud_type,:).Centroid;
            
            if strcmp(data_meta.Sensor,'S_MSI')
                VZAxy = pi/180*mean(single(VZA(orin_cid_all))/100);
                VAAxy = pi/180*mean(single(VAA(orin_cid_all))/100);
            end
            
            % height estimated by neighboring clouds.
            if record_base_h_num>=num_near
                % the centers of already matched clouds
                % current cloud's center
                % the nearest cloud among all matched clouds.
                
                % remove the self cloud heigh
                [nearest_cloud_centers,nearest_dis]=knnsearch(matched_clds_centroid,center_cur,'k',num_near, 'distance','cityblock');% less time-consuming method chebychev
                if cloud_type_cur>num % remove its self for the first 14 clouds when coming back.
                    nearest_cloud_centers(nearest_dis==0)=[];
                end
                
                % get all matched clouds' height.
                record_base_h_tmp=height_clds_recorder(nearest_cloud_centers);
                record_base_h_near=prctile(record_base_h_tmp,100*h_pt);
                h_pre_std=std(record_base_h_tmp);
                clear record_base_h_tmp;
                if h_pre_std>=1000||record_base_h_near <= Min_cl_height||record_base_h_near >= Max_cl_height
                    record_base_h_near=0;
                end
                clear h_pre_std;
            end
           
            dist_pre=0;
            dist_passed=false;
            dist_first=true;
            % all pixels of projected cloud object 
            if numel(orin_cid_all(:))==0
                dist_passed=true;
            else
                cpc_i=center_cur(2);
                cpc_j=center_cur(1);
            end
            
            % indicates the number of the matched cloud shadows for this
            % current cloud.
            num_matched=0;
            
            for base_h=Min_cl_height:step_interval:Max_cl_height % iterate in height (m)
                % Get the true postion of the cloud
                % calculate cloud DEM with initial base height
                if strcmp(data_meta.Sensor,'S_MSI')
                    h=base_h; % have no temperature data. cannot serve as 3D object.
                elseif strcmp(data_meta.Sensor,'L_OLI_TIRS')||...
                        strcmp(data_meta.Sensor,'L_ETM_PLUS')||...
                        strcmp(data_meta.Sensor,'L_TM')
                    h=double(10*(t_obj-temp_obj)/rate_elapse+base_h);% relative to the reference plane. Cloud top's height.
                end
                
                % the height for the bias of the real cloud location should exclude the
                % surface elevation below the cloud object.
                h_bias=h-base_heigh_cloud;% hc-Ec the height between cloud object and its surface.
                if strcmp(data_meta.Sensor,'S_MSI')
                    [tmp_xys(:,1),tmp_xys(:,2)]= getRealCloudPositionS2(orin_xys(:,1),...
                      orin_xys(:,2),h_bias,VZAxy,VAAxy,data_meta.Resolution);
                elseif strcmp(data_meta.Sensor,'L_OLI_TIRS')||...
                        strcmp(data_meta.Sensor,'L_ETM_PLUS')||...
                        strcmp(data_meta.Sensor,'L_TM')
                    sensor_heigh_bias=base_heigh_cloud+dem_base_heigh; % used to exclude the elevation of cloud' surface.
                    [tmp_xys(:,1),tmp_xys(:,2)]=getRealCloudPosition(orin_xys(:,1),...
                        orin_xys(:,2),h_bias,A,B,C,omiga_par,omiga_per,sensor_heigh_bias);
                else
                    error('Only Landsats 4-7, Landsat 8 and Sentinel 2 data can be supported./n');
                end
         
                % shadow moved distance (pixel) to calculate the cloud
                % shadow locaiton.
                % real cloud height relative to the plane.
                i_xy=h/(sub_size*tan(sun_elevation_rad)); 
                XY_type(:,2)=round(tmp_xys(:,1)-i_xy*cos(sun_tazi_rad)); % X is for j,2
                XY_type(:,1)=round(tmp_xys(:,2)-i_xy*sin(sun_tazi_rad)); % Y is for i,1

                clear i_xy;
                % this location is relative to reference plane.
                tmp_j_plane=XY_type(:,2);% col
                tmp_i_plane=XY_type(:,1);% row
                clear XY_type;
                % back project
%                 dim_expd=1000;% 1000 pixels buffer
                % some projected pixels out of observations.
                tmp_i_plane_expd_tmp=tmp_i_plane+dim_expd;
                tmp_j_plane_expd_tmp=tmp_j_plane+dim_expd;
                
                avail_pixels=find(tmp_i_plane_expd_tmp>0&tmp_j_plane_expd_tmp>0&...
                    tmp_i_plane_expd_tmp<=dim_expand(1)&tmp_j_plane_expd_tmp<=dim_expand(2));
                
                tmp_i_plane_expd=tmp_i_plane_expd_tmp(avail_pixels);
                tmp_j_plane_expd=tmp_j_plane_expd_tmp(avail_pixels);
                clear tmp_i_plane_expd_tmp tmp_j_plane_expd_tmp avail_pixels;
                tmp_id_plane_expd=sub2ind(dim_expand,tmp_i_plane_expd,tmp_j_plane_expd); % matched shadow locations
                clear tmp_i_plane_expd tmp_j_plane_expd;
                
                % search the responding locations in real surface (derived 
                % from the relation-lookup table).
                tmp_i_obs=recorderRow(tmp_id_plane_expd);
                tmp_j_obs=recorderCol(tmp_id_plane_expd);
                clear tmp_id_plane_expd tmp_id_plane_expd;
                
                % cloud shadow must be beyond the location of the orgianl cloud.
                if ~dist_passed
                    % the center of cloud shadow in real image.
                    sum_cpmp_i=sum(tmp_i_obs(:));
                    sum_cpmp_j=sum(tmp_j_obs(:));
                    area_cpmp=numel(tmp_j_obs(:));
                    ctmp_i=sum_cpmp_i/area_cpmp;
                    ctmp_j=sum_cpmp_j/area_cpmp;
                    clear sum_cpmp_i sum_cpmp_j area_cpmp;
                    
                    % distance between orginal cloud and its cloud shadow,
                    % Note we ignored the mini bias from view angles here.
                    dist_cur = pdist2([ctmp_j,ctmp_i],[cpc_j,cpc_i],'Euclidean');
%                     dist_cur = floor(dist_cur);
                    clear ctmp_j ctmp_i;
                    if dist_first
                        dist_pre = dist_cur;
                        dist_first = false;
                    else
                        % the distance between the center of cloud and 
                        % cloud shadow over plane decreases
                        if dist_pre >= dist_cur || dist_cur<i_xy_min % should move more than 200 meter high.
                            dist_pre = dist_cur;
                            record_thresh = 0;
%                             record_h=0;
                            continue;
                        else
                            dist_passed = true; % can go 
                        end
                    end
                end
                
                % calculate the similarity for the matched cloud shadow.
                % the id that is out of the entire image
                % out-of-scene pixels should be found. the relationship
                % between locations at plane and DEM is lack.
                out_id=(tmp_i_obs<1|tmp_i_obs>win_height|tmp_j_obs<1|tmp_j_obs>win_width);
                out_all=sum(out_id(:));
                

                % exclude the pixels out of the entire image.
                tmp_ii_obs=tmp_i_obs(out_id==0);
                tmp_jj_obs=tmp_j_obs(out_id==0);
                clear out_id tmp_i_obs tmp_j_obs;
                tmp_id=sub2ind(dim,tmp_ii_obs,tmp_jj_obs); % matched shadow locations
                clear tmp_ii_obs tmp_jj_obs;
                
                out_obs=mask(tmp_id)==0;
                
                id_ex_self = segm_cloud(tmp_id)~=cloud_type;
                % 1st rule: out of obervations; 2nd rule: located in
                % potential shadow or other clouds (exclude the self cloud).

                % Special case #1:
                % for the cloud shadow previoudly matched, the new one
                % cannot reach the boundary and other clouds, which easily
                % result in the overestimation of silimarity.
                match_id_unsure = out_obs | ...
                    (id_ex_self&(data_cloud_potential(tmp_id)==1));
                
                match_id_sure = id_ex_self&data_shadow_potential(tmp_id)==1;
                
                % give half weight to the macthed pixels located in outside and other
                % clouds.
                matched_all=sum(match_id_sure(:))+0.5*sum(match_id_unsure(:))+out_all;
                
                total_all=sum(id_ex_self(:))+out_all;
                
                thresh_match=matched_all/total_all;    
                clear match_id total_all;
                
                
                % used to determine whether the iteration continues or not.
                iter_con=true; % continues as default.
                clear id_ex_self;
                
                if num_matched > 0&&... % already have cloud shadow previously
                        (record_base_h_near > 0 && base_h >= record_base_h_near) % or more than the predicted cloud height.
                     iter_con=false;
                end
                clear pt_unsure;

                % check the matched cloud shadow or not?
                % the following rules are used to decide to continue or not.
                if iter_con
                    if (thresh_match >= Tbuffer*record_thresh)&&...
                            (base_h < Max_cl_height-step_interval)&&...
                            (record_thresh < max_similar)
                        if thresh_match > record_thresh
                            % record max similarity and the corresponding cloud base height.
                            record_thresh=thresh_match;
    %                         record_h=h;
                            record_base_h=base_h;
                        end
                        continue;
                    else 
                        if (record_thresh >= Tsimilar)
                            % successfully find a cloud shadow
                            num_matched=num_matched+1; % indicates one more cloud shadow was found out.
                            % only when expected height available. (record_base_h_near>0)
                            if base_h<record_base_h_near
                                % but allow the searching reach to the neighboring clouds' height
                                if thresh_match>=record_thresh||thresh_match>=max_similar
                                    record_thresh=thresh_match;
    %                                 record_h=h;
                                    record_base_h=base_h;
                                end
                                continue; % much reach the predicted cloud height.
                            end
                        else
                            record_thresh=0;
                            continue;
                        end
                    end
                end
                % 1: continues; 0: not continue and get a cloud shadow
                if num_matched<1
                    break;
                end

                if r_obj>num_pix&&...
                        cloud_type_cur<=num % cannot re-add for the first 14 clouds
                    matched_clds_centroid=[matched_clds_centroid;center_cur];
                    match_clds(cloud_type)=1;
                    height_clds_recorder=[height_clds_recorder;record_base_h];
%                             height_clds_recorder(cloud_type)=record_base_h;
                    record_base_h_num=record_base_h_num+1;
                end
                similar_num(cloud_type)=record_thresh;

                if isequal(data_meta.Sensor,'S_MSI')
                    record_h = record_base_h;
                    h_bias=record_h-base_heigh_cloud;% hc-Ec
                    [tmp_xys_all(:,1),tmp_xys_all(:,2)]= getRealCloudPositionS2(orin_xys_all(:,1),...
                    orin_xys_all(:,2),h_bias,VZAxy,VAAxy,data_meta.Resolution);
                elseif isequal(data_meta.Sensor,'L_OLI_TIRS')||...
                        isequal(data_meta.Sensor,'L_ETM_PLUS')||...
                        isequal(data_meta.Sensor,'L_TM')
                    record_h = double(10*(t_obj-temp_obj_all)/rate_elapse+record_base_h);
                    h_bias=record_h-base_heigh_cloud;% hc-Ec
                    sensor_heigh_bias=base_heigh_cloud+dem_base_heigh;
                    [tmp_xys_all(:,1),tmp_xys_all(:,2)]=getRealCloudPosition(orin_xys_all(:,1),...
                    orin_xys_all(:,2),h_bias,A,B,C,omiga_par,omiga_per,sensor_heigh_bias);
                else
                    error('Only Landsats 4-7, Landsat 8 and Sentinel 2 data can be supported./n');
                end
                clear orin_xys_all;

                i_vir=record_h/(sub_size*tan(sun_elevation_rad));
                
                tmp_XY_type_all(:,2)=round(tmp_xys_all(:,1)-i_vir*cos(sun_tazi_rad)); % X is for col j,2
                tmp_XY_type_all(:,1)=round(tmp_xys_all(:,2)-i_vir*sin(sun_tazi_rad)); % Y is for row i,1

                clear tmp_xys_all i_vir;

                tmp_scol_plane=tmp_XY_type_all(:,2);
                tmp_srow_plane=tmp_XY_type_all(:,1);
                clear tmp_XY_type_all;
                tmp_tmp_i_plane_expd=tmp_srow_plane+dim_expd;
                tmp_tmp_j_plane_expd=tmp_scol_plane+dim_expd;
                clear tmp_srow_plane tmp_scol_plane;

                avail_pixels=find(tmp_tmp_i_plane_expd>0&tmp_tmp_j_plane_expd>0&...
                    tmp_tmp_i_plane_expd<dim_expand(1)&tmp_tmp_j_plane_expd<dim_expand(2));
                tmp_tmp_i_plane_expd=tmp_tmp_i_plane_expd(avail_pixels);
                tmp_tmp_j_plane_expd=tmp_tmp_j_plane_expd(avail_pixels);
                clear avail_pixels;

                % matched shadow locations at plane.
                tmp_tmp_id_plane_expd=sub2ind(dim_expand,tmp_tmp_i_plane_expd,tmp_tmp_j_plane_expd); 
                clear tmp_tmp_i_plane_expd tmp_tmp_j_plane_expd;
                % matched shadow locations at real image (DEM surface).
                tmp_srow=recorderRow(tmp_tmp_id_plane_expd);
                tmp_scol=recorderCol(tmp_tmp_id_plane_expd);
                clear tmp_tmp_id_plane_expd;

                % remove the pixels out of box.
                out_ids=tmp_srow<1|tmp_scol<1|...
                    tmp_srow>win_height|tmp_scol>win_width;
                tmp_srow(out_ids)=[];
                tmp_scol(out_ids)=[];
                clear out_ids;

                tmp_sid=sub2ind(dim,tmp_srow,tmp_scol);
                clear tmp_srow tmp_scol;

                % give shadow_cal=1
%                     data_shadow_matched(tmp_sid)=1;

                if cloud_type_cur>num % re-visit the first 14 clouds.
                    % remove the matched before.
                    data_shadow_matched((data_shadow_matched==cloud_type))=0;
                    % and give new cloud shadow to this.
                end
                data_shadow_matched(tmp_sid)=cloud_type;
                clear tmp_sid;
                clear center_cur;
                break;
            end

        end
        data_cloud_matched=data_cloud_potential;
        data_shadow_matched=uint8(data_shadow_matched>0);
        % remove the cloud.
        data_shadow_matched(data_cloud_matched==1)=0;
    end
end
