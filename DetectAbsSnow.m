function abs_snow = DetectAbsSnow(band_green,band_green_statu,ndsi,psnow,resolution)
%DETECTABSSNOW Select absolute snow/ice pixels using spectral-contextual
%     features.
%
% Syntax
%
%     abs_snow =
%     DetectAbsSnow(band_green,band_green_statu,ndsi,psnow,resolution)
%
% Description
%
%     A spectral-contextual snow index (SCSI) is used to select 100%
%     snow/ice pixels. SCSI10 is computed here, of which 10 indicates 10
%     kilometers-by-10 kilometers window size, that 333-by-333 pixels for
%     Landsat 30 meters image and 501-by-501 pixels for Sentinel-2 20
%     meters image (Fmask runs at this 20 meters). That window size is
%     large enough to capture the various contexts for clouds, but still
%     give us soomth variations for pure snow/ice.
%
% Input arguments
%
%     band_green           Green band.
%     band_green_statu     Statured pixels at green band.
%     psnow                Potential snow/ice.
%     resolution           Spatial resolution of input image.
%
% Output arguments
%
%     abs_snow             Absolute snow/ice.
%
% Example
%
%     abs_snow =
%     DetectAbsSnow(band_green,band_green_statu,ndsi,psnow,resolution)
%
%        
% Author:  Shi Qiu (shi.qiu@uconn.edu)
% Date: 21. January, 2018

    % which is equal to 51*51 % when have snow/ice pixels.
    if sum(psnow(:))> 110889
        % sometimes Nan will appear in this green band.
        band_green_filled = fillmissing(band_green,'constant',0);
        clear band_green;
        % window size is 10 km.
        radius = fix(5000/resolution);
        clear resolution;
        win_size=2*radius+1;
        % compute SCSI.
        % clip to small arrays (16 smalls);
        num_clips = 5;  % 5 cippes
        [size_row, size_colum] =size(band_green_filled);
        num_per_row = round((size_row)/num_clips);
        
        row_nums = zeros([1,num_clips],'double');
        row_nums(1:num_clips-1)=num_per_row;
        % the last one
        row_nums(num_clips) = size_row - num_per_row*(num_clips-1);
        
        scsi = zeros([size_row, size_colum],'double');
        % the first clip
        
        
        for i_row = 1:num_clips
            row_ed = sum(row_nums(1:i_row));
            row_st = row_ed - row_nums(i_row) + 1;
            
            row_st_exp = row_st - radius;
            row_ed_exp = row_ed + radius;
            row_st_small = radius + 1;
            row_ed_small = radius + row_nums(i_row);
            % the first clip
            if row_st_exp < 1
                row_st_exp =1;
                row_st_small =1;
                row_ed_small = row_nums(i_row);% no the first exp
            end
            % the last clip
            if row_ed_exp > size_row
                row_ed_exp = size_row;
%                 row_ed_small = radius + row_nums(i_row);
            end
            scsi_small= stdfilt(band_green_filled(row_st_exp:row_ed_exp,:),...
                     true(win_size)).*(1-ndsi(row_st_exp:row_ed_exp,:));
            scsi(row_st:row_ed,:) = scsi_small(row_st_small:row_ed_small,:);  
            clear scsi_small;
        end
        % scsi=stdfilt(band_green_filled,true(win_size)).*(1-ndsi);
        clear band_green_filled ndsi win_size;
        % only get snow/ice pixels from all potential snow/ice pixels, and 
        % do not select the saturated pixels which may be cloud!
        abs_snow=scsi<9&psnow&(~band_green_statu);
        clear scsi psnow band_green_statu;
    else
        abs_snow=[];
    end
end
