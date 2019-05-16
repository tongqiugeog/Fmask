function cirrus_noml = NormalizaCirrusDEM( mask, idplcd, cirrus, dem )
%NORMALIZACIRRUSDEM normalize Cirrus band using most dark object at Cirrus
%band based on DEMs.
    cirrus_noml=zeros(size(cirrus),'double');
    idclr=idplcd==false&mask; % clear sky pixels
    l_pt=2;
    if isempty(dem) % when no DEMs, we also adjust the cirrus TOA to be around 0;
        cirrus_noml(mask)=cirrus(mask)-prctile(cirrus(idclr),l_pt);
        clear mask cirrus idclr l_pt;
    else
        dem_start=prctile(dem(mask),0.001); dem_end=prctile(dem(mask),99.999); % further exclude errors in DEM.
        clear mask;
        dem_step=100;
        dem_intvl=dem_start:dem_step:dem_end;
        clear dem_start dem_end;
        n_slice=length(dem_intvl);
        dem_LUT=dem_intvl;
        clear dem_intvl;
        cirrus_lowest=0;
        cirrus_lowest_al=0;
        for i=1:n_slice
            % the dem intervals.
            if i==n_slice
                ids_inter=dem>dem_LUT(i)-dem_step/2;
            else
                if i==1
                    ids_inter=dem<dem_LUT(i)+dem_step/2;
                else
                    ids_inter=(dem>dem_LUT(i)-dem_step/2)&(dem<dem_LUT(i)+dem_step/2);
                end
            end
            ids_inter_clr=ids_inter&idclr;
            % we think more than 50 pixels show meaningfull calculation.
            % at the same time, cirrus TOA will be substracted from the
            % previous lowest cirrus of clear pixels.
            if sum(ids_inter_clr(:))>0
                cirrus_lowest=prctile(cirrus(ids_inter_clr),l_pt);
                if cirrus_lowest_al==0 % the first local lowest cirrus value will be given
                    cirrus_lowest_al=cirrus_lowest;
                end
            end
            clear ids_inter_clr;
            
            cirrus_noml(ids_inter)=cirrus(ids_inter)-cirrus_lowest;
            clear ids_inter;
        end
        clear cirrus cirrus_lowest idclr dem_LUT dem_step n_slice;
% %          % when no DEMs, we also adjust the cirrus TOA to be around 0,
% %          % based on the full lowest cirrus value.
% %         cirrus_noml(~dem_have)=cirrus(~dem_have)-cirrus_lowest_al;
    end
    % the normalized cirrus value will be set as 0 when it is negative.
    cirrus_noml(cirrus_noml<0)=0;
end