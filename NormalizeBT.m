function Tempd = NormalizeBT( dim, dem,mask,Temp,ratelaspe_cl,l_pt,h_pt,resl)
%NORMALIZEBT Normalize BT along elevation (DEM) by using a linear model.
%
% Syntax
%
%     Tempd = NormalizeBT( dim, dem,mask,Temp,ratelaspe_cl,l_pt,h_pt )
%
% Description
%
%     A linear model is used to normalize BT along with DEM (Qiu et al.,
%     2017).
%     History:
%     1. Create this function. (1. January, 2017)
%     2. This stratied sampling method sometimes results in no enough
%     samples. If the stratied samples are less than 40,000 (not 50,000), 
%     the stratied sampling method will not be used anymore. (8. March,
%     2018)
%     3. no normalization if no dem data. (20. March, 2018)
%
% Input arguments
%
%     dim            Dim for data.
%     dem            Digital Elevation Model (DEM).
%     Temp           Temperature (BT).
%     ratelaspe_cl   Clear sky (land) pixels, which are used for this
%                    normalization.
%     l_pt           Low level (17.5 percentile).
%     h_pt           High level (81.5 percentile).
%     resl           Spatial resolution (Landsat 30 meters; Sentinel-2 20 meters).
%
% Output arguments
%
%     Tempd          Nomalized Temperature (BT).
%
%        
% Author:  Shi Qiu (shi.qiu@ttu.edu)
% Date: 8. March, 2018
    
    if isempty (dem)
        Tempd = Temp;
    else
        mask(isnan(dem))=0;% exclude nan dem pixel.
        dem_b=double(prctile(dem(mask),0.0001));
        temp_cl=Temp(ratelaspe_cl);
    %     clear ratelaspe_cl;
        temp_min=prctile(temp_cl,l_pt*100);
        temp_max=prctile(temp_cl,h_pt*100);
        clear temp_cl l_pt h_pt;
        cls=(Temp>temp_min&Temp<temp_max)&ratelaspe_cl;
        clear temp_min temp_max ratelaspe_cl;
    %     cls=(temp_cl>temp_min&temp_cl<temp_max);
        data_bt_c_clear=double(Temp(cls));
        data_dem_clear=double(dem(cls));
        clear cls;
        % stratified random samples select
    %     total_sample=50000;
        total_sample=40000;
        ele_strata=300;% meters
        samp_distance=450;% meters
        dem_t=double(prctile(dem(mask),99.999)); % further exclude non dem pixels.
        clear mask temp_max temp_min;
    %         binScatterPlot(data_dem_clear,data_bt_c_clear)
        samples_ids=stratiedSampleHanlder(data_dem_clear,dem_b,dem_t,dim,total_sample,ele_strata,samp_distance,resl);
        clear total_sample ele_strata samp_distance dem_t dim resl;

        data_dem_clear_tmp=data_dem_clear;
        data_bt_c_clear_tmp=data_bt_c_clear;
        clear data_dem_clear data_bt_c_clear;

        data_dem_clear=data_dem_clear_tmp(samples_ids,:);
        data_bt_c_clear=data_bt_c_clear_tmp(samples_ids,:);
        clear data_dem_clear_tmp data_bt_c_clear_tmp samples_ids;

       %% regress
        alpha=0.05;
    %     [b,bint,r,rint,stats]=regress(data_bt_c_clear,[ones(size(data_dem_clear)),data_dem_clear],alpha);
        [b,~,~,~,stats]=regress(data_bt_c_clear,[ones(size(data_dem_clear)),data_dem_clear],alpha);
        clear data_bt_c_clear data_dem_clear;
        rate_lapse=0.0;
        if stats(3)<alpha && double(b(2)) < 0
            rate_lapse=double(b(2));% -0.00
        end

        clear alpha b stats;
        Temp=double(Temp);
        if rate_lapse==0.0
            Tempd=Temp;
        else
            Tempd=Temp-rate_lapse.*double(dem-dem_b);
            clear rate_lapse dem_b;
            Tempd(isnan(dem))=Temp(isnan(dem));% back to orignal values for the no dem pixels.
        end
        clear Temp dem;
    end
end

