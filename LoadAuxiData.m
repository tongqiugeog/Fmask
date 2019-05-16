function [dem_out,slope_out,aspect_out,water_occur] = LoadAuxiData(doc_path,fname,bbox,trgt,saveas_local,varargin)
%LOADAXUIDATA Read DEM, GSWO, slope, and aspect. 
%
% Syntax
%
%     [dem_out,slope_out,aspect_out,water_occur] = 
%     LoadAuxiData(doc_path,fname,bbox,trgt,saveas_local)
%
% Description
%
%     Read DEM and GSWO, and calculate Slope and Aspect data.
%     It will spend about 1 min to load, mosaic DEM and GSWO, and calculate
%     Slope and Aspect. If users use the online DEM, the time will depend
%     on the Internet speed.
%     Create a chache folder at scene directory, that will be affected by
%     each others (especailly in parelel)
%     User can define a DEM path. (by Shi. Jun 14, 2018)
%
% Input arguments
%
%     doc_path     Path of images.
%                  Landsat: the directory where you save the Landsat scene.
%                  Sentinel 2: the directory reaching to
%                  '~/S2*/GRANULE/S2*/IMG_DATA/'.
%     fname        Image name.
%     bbox         Boundary of image. [north, south, west, east]
%     trgt         The GRIDobj used as a targart, which is useful when 
%                  projecting and mosaicing the auxi data.
%     saveas_local Locally save as the auxi data or not.
%                  true: saveas false:do not saveas.
%
% Output arguments
%
%     dem_out      DEM responding to the scene.
%     slope_out    Slope, which will be generated if unvailable.
%     aspect_out   Aspect, which will be generated if unvailable. 
%     water_occur  GSWO responding to the scene.
%
%
%        
% Author:  Shi Qiu (shi.qiu@ttu.edu)
% Date: 22. September, 2017 


    % instant paras
%     dem_path=strcat(pwd,'/AuxiData/GTOPO3010degZIP/'); % cliped by each
%     10 degrees -by- 10 degrees.
    p = inputParser;
    p.FunctionName = 'AuxiParameters';
    addParameter(p,'parallel',0);
    addParameter(p,'userdem','');
     % request user's input
    parse(p,varargin{:});
    isParallel=p.Results.parallel;
    
    userdem = p.Results.userdem; % dem path

    if isdeployed % Stand-alone mode.
%         m_path = ctfroot;
        if ismac||isunix
            % Code to run on Mac plaform
            % Code to run on Linux plaform
            [~, result] = system('echo $PATH');
            m_path = char(regexpi(result, '(.*?):', 'tokens', 'once'));  
        elseif ispc
            % Code to run on PC plaform
            [~, result] = system('path');
            m_path = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
            clear result;
        end
    else % MATLAB mode.
        m_path=mfilename('fullpath');
        m_path(end-12:end)=[]; % remove the '/LoadAuxiData'
    end
    water_path=fullfile(m_path,'AuxiData','GSWO150ZIP');
    dem_path=fullfile(m_path,'AuxiData','GTOPO30ZIP');
    clear m_path;

%     water_path=fullfile('.','AuxiData','GSWO150ZIP');
%     dem_path=fullfile('.','AuxiData','GTOPO30ZIP');
%     cache_path=fullfile(m_path,'Cache');
    cache_path=fullfile(doc_path,'Cache');
    
      % see the Cache is available
    if ~exist(cache_path,'dir')
        %have no the directory, create it here.
        mkdir(cache_path);
%         % cannot make it, give the results here directly.
%         if ~status_mkdir
%             cache_path = doc_path;
%         end
    end
    
    % online

    %% DEM
    % Step 1. Read DEM data at the directory of scene.
    % 	filename_dem=dir([doc_path,fname,'_DEM.tif']);
    % userdem
    dem_gridobj=[];
    if ~isempty(userdem)&&~isempty(trgt)
        % read user's DEM
        userdem_gridobj = GRIDobj(userdem);
        % resample to same extent.
        dem_gridobj = reproject2utm(userdem_gridobj,trgt);
        clear userdem_gridobj userdem;
    end
    % local _DEM.tif
    if isempty (dem_gridobj)
        filename_dem=dir(fullfile(doc_path,'L*_DEM.tif')); % upper 
        if ~isempty(filename_dem)
            dem_gridobj = GRIDobj(fullfile(doc_path,filename_dem(1).name));
        else
            filename_dem=dir(fullfile(doc_path,'L*_dem.tif')); % lower
            if ~isempty(filename_dem)
                dem_gridobj = GRIDobj(fullfile(doc_path,filename_dem(1).name));
            else
                filename_dem=dir(fullfile(doc_path,'S*_DEM.tif')); % upper 
                if ~isempty(filename_dem)
                    dem_gridobj = GRIDobj(fullfile(doc_path,filename_dem(1).name));
                else
                    filename_dem=dir(fullfile(doc_path,'S*_dem.tif')); % lower
                    if ~isempty(filename_dem)
                        dem_gridobj = GRIDobj(fullfile(doc_path,filename_dem(1).name));
                    else
                        dem_gridobj=[];
                    end
                end
            end
        end
    end
%     dem_gridobj=[];
    clear filename_dem;
    % Step 2. Locally load GTOPO30 DEM (1km) or Download SRTM DEM data according to users' requirements.
    if isempty (dem_gridobj)&&~isempty(trgt)
        islocal=true;
        if islocal
            dem_data=MosaicLocalDEM(dem_path,cache_path,trgt,bbox);
%             dem_data = MosaicLocalAuxiData('DEM',dem_path,trgt,bbox);
            if ~isempty(dem_data)&&sum(sum(dem_data(:)))>0 % all 0 will be empty.
                % sometimes there may be over water totally.
                dem_gridobj = GRIDobj(trgt);
                dem_gridobj.name=[fname,'_DEM'];
                dem_gridobj.Z = dem_data;
            
%             tic
%             while sum(sum(isnan(dem_gridobj.Z )))>0
%                 dem_gridobj = inpaintnans(dem_gridobj,'fill');% imfill nan
%             end
%             toc
                if saveas_local
                    GRIDobj2geotiff(dem_gridobj,fullfile(doc_path,[dem_gridobj.name,'.tif']));
                end
%                 fprintf('GTOPO30 (1km DEM), ');
            else
%                 fprintf('non-DEM, ');
            end
            clear dem_data;
%             fprintf('Local GTOPO30 (1km DEM) will be used.\n');
        else
            demtype = 'SRTMGL3';
            dem_gridobj = LoadOnlineDEM(fname,demtype,bbox,trgt);
            if ~isempty(dem_gridobj)
                if saveas_local
                    dem_gridobj.name=[fname,'_DEM'];
                    GRIDobj2geotiff(dem_gridobj,fullfile(doc_path,[dem_gridobj.name,'.tif']));
                end
%                 fprintf('Online DEM (',demtype,'),\n');
            else
%                 fprintf('non-DEM (',demtype,'), ');
            end
            
        end
    elseif ~isempty (dem_gridobj)
%         fprintf('DEM in the scene directory,\n');
%         fprintf('The DEM in the directory of scene will be used.\n');
    end
    
    if ~isempty(dem_gridobj)
        dem_out=single(dem_gridobj.Z);
    else
        dem_out=[];
    end
    
    % Step 3. Calculate Slope.
    slope_out = LoadAuxiItem(doc_path,fname,'SLOPE');
    if isempty (slope_out)
        if ~isempty (dem_gridobj)
            slope_gridobj=arcslope(dem_gridobj,'deg');
            slope_out=single(slope_gridobj.Z);
            if saveas_local
                slope_gridobj.name=[fname,'_SLOPE'];
                GRIDobj2geotiff(slope_gridobj,fullfile(doc_path,[slope_gridobj.name,'.tif']));
            end
            clear slope_gridobj;
%             fprintf('SLOPE, ');
%             clear slope_gridobj;
        else
%             fprintf('non-SLOPE, ');
        end
    else
%         fprintf('SLOPE in the scene directory,\n');
    end
    % Step 4. Calculate Aspect.
    aspect_out = LoadAuxiItem(doc_path,fname,'ASPECT');
    if isempty (aspect_out)
        if ~isempty (dem_gridobj)
            aspect_gridobj = aspect(dem_gridobj);
            aspect_out = single(aspect_gridobj.Z);
            if saveas_local
                aspect_gridobj.name=[fname,'_ASPECT'];
                GRIDobj2geotiff(aspect_gridobj,fullfile(doc_path,[aspect_gridobj.name,'.tif']));
            end
%             fprintf('ASPECT, ');
            clear aspect_gridobj;
        else
%             fprintf('non-ASPECT, ');
        end
    else
%         fprintf('ASPECT in the scene directory,\n');
    end
    clear dem_gridobj;
    
    %% WATER
    water_occur = LoadAuxiItem(doc_path,fname,'WATER');

    if isempty (water_occur)&& ~isempty(trgt)
        water_occur = MosaicLocalAuxiData('WATER',water_path,cache_path,trgt,bbox);
        if ~isempty(water_occur)
%             fprintf('and GSWO (150m water map).\n');
%             fprintf('Local GSWO (150m water map) will be used.\n');
            if saveas_local
                wt_gridobj = GRIDobj(trgt);
                wt_gridobj.name=[fname,'_WATER'];
                wt_gridobj.Z = water_occur;
                GRIDobj2geotiff(wt_gridobj,fullfile(doc_path,[wt_gridobj.name,'.tif']));
                clear wt_gridobj;
            end
        else
%             fprintf('and non-GSWO.\n');
        end
    else
%         fprintf('and GSWO in the scene directory.\n');
    end
    clear trgt;
    if ~isempty (water_occur)
        water_occur(water_occur==255)=100;% 255 is 100% ocean.
    end
    if ~isParallel
        % remove the cache folder in the scene directory.
        if exist(cache_path,'dir')
            rmdir(cache_path,'s');
        end
    end
end


function auxi_data = LoadAuxiItem(doc_path,fname,name_match)
% LoadAuxiItem read local auxi data first.

% 	filename=dir([doc_path,'x*_',name_match,'.tif']);
	filename=dir(fullfile(doc_path,[fname,'_',name_match,'.tif']));
    if ~isempty(filename)
        auxi_data=imread(fullfile(doc_path,filename(1).name));
%         fprintf('Successfully load %s at %s.\n',name_match,fname);
    else
        auxi_data=[];
%         fprintf('Fail to locally find out %s at %s.\n',name_match,fname);
    end
end

function dem = LoadOnlineDEM(fname,demtype,bbox,trgt)
    try
%         fprintf('DEM is downloading.\n');
%         cache_path=fullfile(pwd,'Cache');
        tempfile_path=fullfile(cache_path,[fname, '_DEM.tif']);
        dem_temp=Readopentopo('filename',tempfile_path,...
            'north',bbox(1),...
            'south',bbox(2),...
            'west',bbox(3),...
            'east',bbox(4),...
            'demtype',demtype,...
            'deletefile',true);
        if ~isequal(dem_temp.size,trgt.size) % when DEM is NOT same as Landsat, CLIP it.
            dem = reproject2utm(dem_temp,trgt,'method','nearest');
        else
            dem = dem_temp;
        end
%         fprintf('DEM was successfully downloaded at %s.\n',fname);
    catch
        dem=[];
%         fprintf('DEM was unsuccessfully downloaded at %s.\n',fname);
    end
end

function auxi_data=MosaicLocalDEM(auxi_path,cache_path,target_gridobj,bbox)
    latlim=[bbox(2) bbox(1)];
    lonlim=[bbox(3) bbox(4)];
    clear bbox;
%     bbox=[north,south,west,east];
%     [southern_limit northern_limit]
%     [western_limit eastern_limit]

    tileNames = gtopo30s(latlim,lonlim);
    clear latlim  lonlim;
    if isempty(tileNames) % no tiles, return empty.
        auxi_data=[]; % set empty layer.
        return;
    end
%     auxi_data = zeros(target_gridobj.size,'double');  % dem mask
    auxi_data = NaN(target_gridobj.size,'double');  % dem mask
%     mask_dem = createmask(target_gridobj);

    for i=1:length(tileNames)
        tile_name=['gt30',lower(char(tileNames(i)))];
        if ~isempty(tile_name)
            % unzip files.
            filename = fullfile(auxi_path,[tile_name,'.zip']);
            if exist(filename,'file')
                unzip(filename,cache_path);
                filename = fullfile(cache_path,[tile_name,'.tif']);
                % reproject to UTM same with target image.
                auxi_data_part_tmp = GRIDobj(filename);
                clear filename;

                % In the gtopo30 DEM, ocean areas have been masked as "no data"
%                 ids_nan=auxi_data_part_tmp.Z==-9999;
%                 auxi_data_part_tmp.Z(ids_nan)=0;% give 0 to ocean.
                % further give 0 to nan.
                auxi_data_part_tmp.Z=fillmissing(auxi_data_part_tmp.Z,'constant',0);

% %                 % recontruct the refmat to 3-by-2 matrix because of the orginal RGB.
% %                 LON11=min(auxi_data_part_tmp.georef.SpatialRef.LongitudeLimits);
% %                 LAT11=max(auxi_data_part_tmp.georef.SpatialRef.LatitudeLimits);
% %                 % Convert to a geographic raster reference object:
% %                 DLON=auxi_data_part_tmp.cellsize;
% %                 DLAT=0-auxi_data_part_tmp.cellsize;
% %                 auxi_data_part_tmp.refmat = makerefmat(LON11, LAT11, DLON, DLAT);
                auxi_data_tmp = reproject2utm(auxi_data_part_tmp,target_gridobj);
                clear auxi_data_part_tmp;
                auxi_data=max(auxi_data,auxi_data_tmp.Z);
                clear auxi_data_tmp;
% %                 delete(filename); % delete the files
            end
            clear tile_name;
        end
    end
    clear tileNames;
    % if using more than 2 tiles, the mosaic method maybe result in some
    % Nan values at lines. fast fill then.
    if sum(sum(isnan(auxi_data)))>0
        % see along-row or along-col
        [row_nan,col_nan]=find(isnan(auxi_data));
        row_diff=max(row_nan)-min(row_nan);
        col_diff=max(col_nan)-min(col_nan);
        % along row and col.
        if row_diff>col_diff
            auxi_data=fillmissing(auxi_data,'linear',2);% by row.
        else
            auxi_data=fillmissing(auxi_data,'linear',1);% by column.
        end
        % https://www.mathworks.com/help/matlab/ref/fillmissing.html?searchHighlight=fillmissing&s_tid=doc_srchtitle
    end
    % when have nan values in DEM,this will be ocean. all as 0.
 	auxi_data=fillmissing(auxi_data,'constant',0);
end

function auxi_data = MosaicLocalAuxiData(auxi_type,auxi_path,cache_path,target_gridobj,bbox)

    switch auxi_type
       case 'DEM'
            pname_1th='gt30';
       case 'WATER'
            pname_1th='occurrence';
       otherwise
            error('The input of auxiliary type is incorrect, which should be ''DEM'' or ''WATER''.');
    end
    
    % target_image
    % get the up-left cornner's name for water title
    [tiles_used, ul_name, ur_name, ll_name, lr_name] = ConvertBBox2TileName(bbox,pname_1th);
    clear bbox pname_1th;
    if sum(tiles_used(:))==0 % no tiles, return empty.
        auxi_data=[]; % set empty water layer.
        return;
    end
    
    auxi_data = zeros(target_gridobj.size,'double');  % water mask
    % loop all water titles covering the Landsat or Sentinel image.
    for i=1:4
        tile_name='';
        switch i
            case 1
               if tiles_used(i)
                    tile_name=ul_name;
               end 
            case 2
               if tiles_used(i)
                    tile_name=ur_name;
               end 
            case 3
               if tiles_used(i)
                    tile_name=ll_name;
               end 
            case 4
               if tiles_used(i)
                    tile_name=lr_name;
               end
        end
        if ~isempty(tile_name)
            % unzip files.
            filename = fullfile(auxi_path,[tile_name,'.zip']);
%             cache_path=fullfile(pwd,'Cache');
            if exist(filename,'file')
                unzip(filename,cache_path);
                filename = fullfile(cache_path,[tile_name,'.tif']);
%             filename = [water_path,tile_name,'.tif'];
            % reproject to UTM same with target image.
                auxi_data_part_tmp = GRIDobj(filename);
                clear filename;
                % recontruct the refmat to 3-by-2 matrix because of the orginal RGB.

%                 LON11=min(auxi_data_part_tmp.georef.SpatialRef.LongitudeLimits);
%                 LAT11=max(auxi_data_part_tmp.georef.SpatialRef.LatitudeLimits);
%                 % Convert to a geographic raster reference object:
% %                 rasterSize = auxi_data_part_tmp.size;
%                 DLON=auxi_data_part_tmp.cellsize;
%                 DLAT=0-auxi_data_part_tmp.cellsize;
% %                 cell_resol=10/rasterSize(1);
% %                 DLON=cell_resol;
% %                 DLAT=0-cell_resol;
%                 auxi_data_part_tmp.refmat = makerefmat(LON11, LAT11, DLON, DLAT);

                auxi_data_part_tmp.Z=fillmissing(auxi_data_part_tmp.Z,'constant',255);% 255 or No data indicates ocean area.
%                 auxi_data_part_tmp.Z(isnan(auxi_data_part_tmp.Z))=255;% 255 or No data indicates ocean area.
                auxi_data_tmp = reproject2utm(auxi_data_part_tmp,target_gridobj);
                clear auxi_data_part_tmp;
%                 auxi_data_part_tmp.Z=fillmissing(auxi_data_part_tmp.Z,'constant',0);% no data indicates no sure whether there is water.
                auxi_data=max(auxi_data,auxi_data_tmp.Z);
                clear auxi_data_tmp;
%                 delete(filename); % delete the files
            else
                tiles_used(i)=0;% no data file!
            end
        end
    end
    % see it has or not again
    if sum(tiles_used(:))==0 % no tiles, return empty.
        auxi_data=[]; % set empty.
        return;
    end
    
    if sum(tiles_used(:))>1% masic will lead to none data. 
        auxi_data = medfilt2(auxi_data); % filter it.
    end
        
    clear tiles_used;
    
    switch auxi_type
       case 'DEM'
        auxi_data=int16(auxi_data); % convert to intergr with 16 bit for DEM, can be negative and positive
       case 'WATER'
        auxi_data=uint8(auxi_data); % convert to intergr with 8 bit for water, only be positive.
    end
    clear auxi_type;
end

function [tiles_used, ul_name, ur_name, ll_name, lr_name] = ConvertBBox2TileName(bbox,pname_1th)
    tiles_used=zeros(1,4);
%    bbox=[north,south,west,east];
% e.g., 21.2993   19.1468   77.7836   80.0166
% 30N 19N 70E 90E
    size_tile=10;% unit: degree.
%     bbox=[21.2993   19.1468   77.7836   80.0166];
    % to negative infinity
    bbox_floor=size_tile.*floor(bbox./size_tile); %e.g., floor: floor(-1.3)=-2; floor(1.3)=1;
    % 20    10    70    80
    % to positive infinity
    bbox_ceil=size_tile.*ceil(bbox./size_tile); %e.g., ceil: ceil(-1.3)=-1; ceil(1.3)=2;
    clear size_tile bbox;
    % 30    20    80    90
    
% %     % four corners for each image.
% %     ul_loc=[bbox(1) bbox(3)];
% %     ur_loc=[bbox(1) bbox(4)];
% %     ll_loc=[bbox(2) bbox(3)];
% %     lr_loc=[bbox(2) bbox(4)];
    % find the Granule's top-left corner.
    ul_loc=[bbox_ceil(1) bbox_floor(3)];
    ur_loc=[bbox_ceil(1) bbox_floor(4)];
    ll_loc=[bbox_ceil(2) bbox_floor(3)];
    lr_loc=[bbox_ceil(2) bbox_floor(4)];
    clear bbox_ceil bbox_floor;
%     pname_1th='occurrence';
    ul_pname_2th=convertNumel2Char(ul_loc);
    ur_pname_2th=convertNumel2Char(ur_loc);
    ll_pname_2th=convertNumel2Char(ll_loc);
    lr_pname_2th=convertNumel2Char(lr_loc);
    clear ul_loc ur_loc ll_loc lr_loc;
    
    
    ul_name=[pname_1th,ul_pname_2th];
    ur_name='';
    ll_name='';
    lr_name='';
    % remove repeat tiles.
    tiles_used(1)=1; % always have 1 or more than 1 tile.
    
    if ~isequal(ul_pname_2th,ur_pname_2th)
        ur_name=[pname_1th,ur_pname_2th];
        tiles_used(2)=1;
    end
    if ~isequal(ul_pname_2th,ll_pname_2th)&&~isequal(ur_pname_2th,ll_pname_2th)
        ll_name=[pname_1th,ll_pname_2th];
        tiles_used(3)=1;
    end
    if ~isequal(ul_pname_2th,lr_pname_2th)&&~isequal(ur_pname_2th,lr_pname_2th)&&~isequal(ll_pname_2th,lr_pname_2th)
        lr_name=[pname_1th,lr_pname_2th];
        tiles_used(4)=1;
    end
end

function pname_2th=convertNumel2Char(corner)
    pname_2th='_';
    if corner(2)>=0
        pname_2th=[pname_2th,num2str(abs(corner(2))),'E_'];
    else
        pname_2th=[pname_2th,num2str(abs(corner(2))),'W_'];
    end
    if corner(1)>=0
        pname_2th=[pname_2th,num2str(abs(corner(1))),'N'];
    else
        pname_2th=[pname_2th,num2str(abs(corner(1))),'S'];
    end
end
