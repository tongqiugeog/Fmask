function samples_ids=stratiedSampleHanlder(data_dem_clear,dem_b,dem_t,dim,total_sample,ele_strata,distance,resolution)
%STRATIEDSAMPLEHANLDER Further select clear sky (land) pixels to normalize
%BT.
%
% Syntax
%
%     samples_ids=
%     stratiedSampleHanlder(data_dem_clear,dem_b,dem_t,dim,total_sample,ele_strata,distance)
%
% Description
%
%     Stratied sampling method is used by Qiu et al., (2017).
%     History:
%     1. Create this function. (1. January, 2017)
%     2. This stratied sampling method sometimes results in no enough
%     samples. (8. March, 2018)
%
% Input arguments
%
%     data_dem_clear     Clear sky pixels' DEM.
%     dem_b              Basic elevation.
%     dem_t              Highest elevation.
%     dim                Dim for data.
%     total_sample       All clear sky samples (40,000).
%     ele_strata         Stratied sampling along elevation (300 meters).
%     distance           Minmum distance among different samples(450 meters).
%                        Distance rule will be removed if distance = 0.
%     resolution         Spatial resolution (Landsat 30 meters; Sentinel-2 20 meters).
% Output arguments
%
%     Tempd          Nomalized Temperature (BT).
%
%        
% Author:  Shi Qiu (shi.qiu@ttu.edu)
% Date: 8. March, 2018


% %     num_clear_sky_pixels=numel(data_dem_clear);
% %     % no enough clear-sky pixels afther removing the pixels out of the upper
% %     % and lower levels.more than 1/4 total samples (10,000) can be used for 
% %     % estimating lapse rate and c in topo correction.
% %     if num_clear_sky_pixels<total_sample/4 
% %         samples_ids=1:num_clear_sky_pixels;
% %         return;
% %     end

    samp_dis= round(distance/resolution); % 450 m(15 30 m pixels)
    clear distance resolution;
    % compute the number of available strata
    strata_avail=[];
    for i_dem=dem_b:ele_strata:dem_t
        dem_clear_index_tmp=data_dem_clear>=i_dem&data_dem_clear<i_dem+ele_strata;
        if sum(dem_clear_index_tmp)>0
%                 num_strata=num_strata+1;
            strata_avail=[strata_avail;1];
        else
            strata_avail=[strata_avail;0];
        end
        clear dem_clear_index_tmp;
    end
    % equal samples in each stratum
    num_sam_each_strata=round(total_sample/sum(strata_avail));
    clear total_sample;
    samples_ids=[];% to restore samples selected.
    % loop each strata and select samples
    loop_i=1;
    for i_dem=dem_b:ele_strata:dem_t
        if strata_avail(loop_i)
%                 dem_clear_index_tmp=data_dem_clear>=i_dem&data_dem_clear<i_dem+ele_strata;
            % find all clear-sky pixels in subgroup.
            samples_ids_tmp=find(data_dem_clear>=i_dem&data_dem_clear<i_dem+ele_strata);
            % randomly selection.
            samples_ids_rand=samples_ids_tmp(randperm(numel(samples_ids_tmp))); 
            clear samples_ids_tmp;
            num_tmp=size(samples_ids_rand,1);
            if samp_dis==0 % no distance rule
                num_max_tmp=num_sam_each_strata;
                if num_tmp>num_max_tmp
                    num_tmp=num_max_tmp;
                end
                clear num_max_tmp;
                samples_ids_rand_tmp=samples_ids_rand(1:num_tmp);
                clear num_tmp;
                % store data
                samples_ids=[samples_ids;samples_ids_rand_tmp];
                clear samples_ids_rand samples_ids_rand_tmp;
            else
                num_max_tmp=num_sam_each_strata*samp_dis;
                if num_tmp>num_max_tmp
                    num_tmp=num_max_tmp;
                end
                clear num_max_tmp;
                samples_ids_rand_tmp=samples_ids_rand(1:num_tmp);
                clear num_tmp;
                samples_ids_rand=samples_ids_rand_tmp;
                clear samples_ids_rand_tmp;
                all_num_tmp=size(samples_ids_rand);
                [i_tmp,j_tmp]=ind2sub(dim,samples_ids_rand);
                ij_tmp=[i_tmp,j_tmp];
                clear i_tmp j_tmp;
                % removing the clear-sky pixels of which distance lt 15. 
                idx_lt15 = rangesearch(ij_tmp,ij_tmp,samp_dis, 'distance','cityblock');
                recorder_tmp=zeros(all_num_tmp,'uint8')+1;
                for i_idx=1:all_num_tmp
                    if recorder_tmp(i_idx)>0 % when this label is available.
                        recorder_tmp(cell2mat(idx_lt15(i_idx)))=0;
                        recorder_tmp(i_idx)=1;% reback the current label as 1.
                    end 
                end
                clear all_num_tmp idx_lt15;
                idx_used=find(recorder_tmp==1);
                clear recorder_tmp;
                num_tmp=size(idx_used,1);
                num_max_tmp=num_sam_each_strata;
                if num_tmp>num_max_tmp
                    num_tmp=num_max_tmp;
                end
                idx_used=idx_used(1:num_tmp);
                clear num_tmp num_max_tmp;
                % store data
                samples_ids=[samples_ids;samples_ids_rand(idx_used)];
                clear samples_ids_rand idx_used;
            end
        end
        loop_i=loop_i+1;
    end
    clear samp_dis num_sam_each_strata;

% %     % no enough clear-sky pixels afther removing the pixels out of the upper
% %     % and lower levels.more than 1/4 total samples (10,000) can be used for 
% %     % estimating lapse rate and c in topo correction.
% %     if numel(samples_ids)<total_sample/4
% %         samples_ids=1:num_clear_sky_pixels;
% %         return;
% %     end
end
       