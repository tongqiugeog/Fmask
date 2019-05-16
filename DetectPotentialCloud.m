function [sum_clr,cloud,idused,t_templ,t_temph]=DetectPotentialCloud(...
      data_meta,mask,water,data_toabt, dem, ndvi,ndsi,ndbi,idplcd,...
      whiteness,HOT,wpt,cldprob)
%DETECTPOTENTIALCLOUD Detect potential clouds using scene-based method.
%
% Syntax
%
%     [sum_clr,cloud,idlnd,t_templ,t_temph]=
%     DetectPotentialCloud(data_meta,mask,water,data_toabt, dem, ndvi,ndsi,
%     ndbi,idplcd,whiteness,HOT,wpt,cldprob)
%
% Description
%
%     Several cloud probabilities are combinated together to capture cloud
%     features of white, bright, cold, and/or high.
%
% Input arguments
%
%     data_meta      Metadata including [row size, column size].
%     mask           Observation mask (outside).
%     water          Water mask.
%     data_toabt     TOA reflectance and BT.
%     dem            DEM data.
%     ndvi           NDVI.
%     ndsi           NDSI.
%     ndbi           NDBI.
%     idplcd         Absolute clear sky pixels.
%     whiteness      Whitness.
%     HOT            Data derived from the HOT transform.
%     wpt            Weight of thin probability (0.3 for Landsat 8 and 
%                    0.5 for Sentinel 2).
%     cldprob        Cloud probability threshold to segment clouds from
%                    surface.
%
% Output arguments
%
%     sum_clr        The total number of clear sky pixels.
%     cloud          Potential clouds.
%     idlnd          Clear sky land pixels.
%     t_templ        Low level temperature (78.5 percentile).
%     t_temph        High level temperature (81.5 percentile).
%
%
%        
% Author:  Shi Qiu (shi.qiu@uconn.edu)
% Date: 20. January, 2018

    % inputs: BandCirrus BandBT BandSWIR1 SatuGreen SatuRed
    cloud = zeros(data_meta.Dim,'uint8');  % cloud mask
    %% Constants Parameters
    l_pt=0.175; % low percent
    h_pt=1-l_pt; % high percent
    
    %% Step 2: calcualte cloud probability
    % select clear sky pixels for land and water, respectively.
    idclr=idplcd==false&mask==1;
    sum_clr=sum(idclr(:));
    idlnd = idclr&water==false;
    idwt = idclr&water==true;%&data(:,:,6)<=300;
    
    t_templ = 0;
    t_temph = 0;
    idused = [];
  % 99.9% TO 99.99%
    if sum_clr <= 40000 % when potential cloud cover less than 0.1%, directly screen all PCPs out.
        cloud(idplcd==true)=1; % all cld
        cloud(mask==0)=0;
    else
        %%%%%%%%%%% thin cloud prob for both water and land%%%%%%%%%%%
        prob_thin = 0; % there is no contribution from the new bands
        if ~isempty(data_toabt.BandCirrus) % Landsat 4~7
            prob_thin = probThin(data_toabt.BandCirrus);
            data_toabt.BandCirrus = [];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%cloud prob over water%%%%%%%%%%%%%%%%%%%
        wprob_temp=1;
        if ~isempty(data_toabt.BandBT)
            wprob_temp = probwTemperature(data_toabt.BandBT,idwt,h_pt);
        end
        
        wprob_brightness = probwBrightness(data_toabt.BandSWIR1);
        data_toabt.BandSWIR1 = [];
        
        %%%%%%%%%%%%%%%%%%%%%%cloud prob over land%%%%%%%%%%%%%%%%%%%%
        lndptm=100*sum(idlnd(:))/sum(mask(:));
        if lndptm >= 0.1
            idused=idlnd;
        else % when having no enough clear land pixels, we used all PCPs to calculate clear land basics.
            idused=idclr;
        end
        clear lndptm;
        
        lprob_temp=1;
        lprob_brightness=1;
        
        
        if ~isempty(data_toabt.BandBT) % if have BT. normalize it using DEMs and use it to calcualte temperature probability.
            data_toabt.BandBT = NormalizeBT( data_meta.Dim,dem,mask,data_toabt.BandBT,idused,l_pt,h_pt,data_meta.Resolution(1));
            [lprob_temp, t_templ,t_temph]=problTemperature(data_toabt.BandBT,idused,l_pt,h_pt);
        else % if have no BT, use HOT probability instead of temperature probability.
            lprob_brightness = problBrightness(HOT,idused,l_pt,h_pt);
            clear HOT l_pt;
        end
%         clear idused;
        
        lprob_vari = problSpectralVaribility(ndvi,ndsi,ndbi,whiteness,data_toabt.SatuGreen,data_toabt.SatuRed);
        clear ndvi ndsi whiteness;

        %%%%%%%%%%%%%%%%%%%%%%%%%%final clouds%%%%%%%%%%%%%%%%%%%%%%%
        % [Final prob mask (water)]
        wprob_final=wprob_temp.*wprob_brightness + wpt.*prob_thin; % cloud over water probability
        clear wprob_temp wprob_brightness;
        wprob_final=100.*wprob_final; % convert percentage
        wclr_h=prctile(wprob_final(idwt),100*h_pt);
        clear idwt;
        % Final prob mask (land)
        lprob_final=lprob_temp.*lprob_vari.*lprob_brightness + wpt.*prob_thin; % cloud over land probability
%         clear lprob_temp lprob_vari lprob_brightness prob_thin wpt;
        clear lprob_temp lprob_vari prob_thin wpt;
        lprob_final=100.*lprob_final; % convert percentage
        clr_h=prctile(lprob_final(idlnd),100*h_pt);
        clear h_pt;
        
        wclr_max=wclr_h+cldprob;% dynamic threshold (water)
        clr_max=clr_h+cldprob;% dynamic threshold (land)
        clear cldprob;
        clear idclr;

        % all potential clouds
% %         cloud(idclr==true)=-1;
       
        id_final_cld = SegementClouds(idplcd,water,data_toabt.BandBT,t_templ,lprob_final,wprob_final,clr_max,wclr_max);
        clear clr_max wclr_max;
        cloud(id_final_cld)=1;
        clear id_final_cld;
        
        cloud(mask==0)=0;
        clear mask;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%final clouds%%%%%%%%%%%%%%%%%%%%%%%
function id_final_cld = SegementClouds(idplcd,water,bt,t_templ,lprob,wprob,clr_max,wclr_max)
    % final clouds
    id_final_cld=idplcd==true&((lprob>clr_max&water==0)|...% cloud over land
            (wprob>wclr_max&water==1));% thin cloud over water
    if ~isempty(bt) % if have BT.
        id_final_cld=id_final_cld|(bt<t_templ-3500);% extremly cold cloud
    end
end

