function [data_meta,data_toabt,angles_view,trgt] = LoadData(path_data,sensor,InputFile,main_meta)
%LOADDATA Read all used bands and metadata for Landsats 4-8 and Sentinel 2.
%
% Syntax
%
%     [data_meta,data_toabt,trgt] = LoadData(path_data)
%
% Description
%
%     For Landsat Level-1 images, Digital Number (DN) values are converted 
%     to TOA reflectances and BT (Celsius degree) by using the LEDAPS 
%     atmosphere correction tool (Masek et al., 2006). This function is 
%     derived from Fmask 3.3 for Landsat.
%     For Sentinel 2 Level-1C images, TOA reflectances are provided with a 
%     spatial resolution of 10, 20 and 60 meters of the different spectral 
%     bands, that will be resampled into 20 meters at Fmask 4.0 rountine.
%     This read function is derived from Fmask 3.3 for Sentinel 2.
%
% History
%
%     1. Create this function. (17. May, 2017) 
%
% Input arguments
%
%     path_data    Path of images.
%                  Landsat: the directory where you save the Landsat scene.
%                  Sentinel 2: the directory reaching to '~/S2*/GRANULE/S2*/IMG_DATA/'
%
% Output arguments
%
%     data_meta    Metadata
%     data_toabt   TOA reflectances (x10000) and BT(x100 Celsius degree) data.
%     trgt         The GRIDobj used as a targart, which is useful when 
%                  projecting and mosaicing the auxi data.
%
% See also: nd2toarbt nd2toarbt_msi
%
%        
% Author:  Shi Qiu (shi.qiu@ttu.edu)
% Date: 17. May, 2017 

    trgt=[];
    angles_view=[];
    if strcmp(sensor , 'S_MSI' )
        
%         fprintf('Load TOA reflectances from the Level-1C product.\n');
        % Load data for Sentinel 2
        [~,data,trgt,dim,bbox,ul,zen,azi,zc,Angles,satu_B1,satu_B2,satu_B3,resolu]=nd2toarbt_msi(InputFile);
        %% Input data class defination
        data_meta=ObjMeta;
        data_toabt=ObjTOABT;
        %% Parameters

        norln=strread(main_meta.name,'%s','delimiter','.'); 
        n_name=char(norln(1));
        data_meta.Name=n_name;
        data_meta.Sensor=sensor;
        data_meta.Dim=dim;
        data_meta.UL=ul;
        data_meta.Zen=zen;
        data_meta.Azi=azi;
        data_meta.ZC=zc;
        data_meta.Resolution=resolu;
        data_meta.BBox=bbox;
        
        fmask_dir='FMASK_DATA';
        % see the FMASK_DATA is available
        if ~exist(fullfile(InputFile.pathh,fmask_dir),'dir')
            %have no the Fmask directory, create it here.
            status_mkdir = mkdir(fullfile(InputFile.pathh,fmask_dir));
            % cannot make it, give the results here directly.
            if ~status_mkdir
                fmask_dir='';
            end
        end
        data_meta.Output=fmask_dir;
        clear fmask_dir;

        %% TOA
        data_toabt.BandBlue=data(:,:,1);
        data_toabt.BandGreen=data(:,:,2);
        data_toabt.BandRed=data(:,:,3);
        data_toabt.BandNIR=data(:,:,4); % band 8A
        data_toabt.BandSWIR1=data(:,:,5);
        data_toabt.BandSWIR2=data(:,:,6);
        data_toabt.BandCirrus=data(:,:,7);
        
        data_toabt.BandVRE3 = data(:,:,8); % band 07
        data_toabt.BandNIR8 = data(:,:,9); % band 08
        
        %% View angles
        angles_view=Angles;

        %% Saturated at visible bands
        data_toabt.SatuBlue=satu_B1;
        data_toabt.SatuGreen=satu_B2;
        data_toabt.SatuRed=satu_B3;
    else
%         fprintf('Calculate TOA reflectances and BT from the Collection 1 product.\n');
        % Load data for Landsat 4-8
        [Temp,data,trgt,dim,ul,bbox,zen,azi,zc,satu_B1,satu_B2,satu_B3,resolu]=nd2toarbt(path_data,main_meta.name);

        %% Input data class defination
        data_meta=ObjMeta;
        data_toabt=ObjTOABT;
        %% Parameters

        % reedit dir_im
        norln=strread(main_meta.name,'%s','delimiter','.'); 
        n_name=char(norln(1));
        data_meta.Name=n_name(1:end-4);
        data_meta.Sensor=sensor;
        data_meta.Dim=dim;
        data_meta.UL=ul;
        data_meta.Zen=zen;
        data_meta.Azi=azi;
        data_meta.ZC=zc;
        data_meta.Resolution=resolu;
        data_meta.BBox=bbox;
        data_meta.Output='';

        %% TOA
        if strcmp(sensor,'L_ETM_PLUS')||strcmp(sensor,'L_TM')
            data_toabt.BandBlue=data(:,:,1);
            data_toabt.BandGreen=data(:,:,2);
            data_toabt.BandRed=data(:,:,3);
            data_toabt.BandNIR=data(:,:,4);
            data_toabt.BandSWIR1=data(:,:,5);
            data_toabt.BandSWIR2=data(:,:,6);
        else
            if strcmp(sensor,'L_OLI_TIRS')
                data_toabt.BandBlue=data(:,:,1);
                data_toabt.BandGreen=data(:,:,2);
                data_toabt.BandRed=data(:,:,3);
                data_toabt.BandNIR=data(:,:,4);
                data_toabt.BandSWIR1=data(:,:,5);
                data_toabt.BandSWIR2=data(:,:,6);
                data_toabt.BandCirrus=data(:,:,7);
            end
        end

        %% BT
        
%         fprintf('BT, ');
        data_toabt.BandBT=Temp;
        %% Saturated at visible bands
        data_toabt.SatuBlue=satu_B1;
        data_toabt.SatuGreen=satu_B2;
        data_toabt.SatuRed=satu_B3;
    end
end

