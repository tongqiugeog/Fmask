function [num_all_images, sensors, paths, info_count_text] = CheckImagesPath(path_data)
%CHECKIMAGEPATH Ensure the input path is right to find a Landsats 4-8, and
%Sentinel 2 image(s).
% path_data - the input path

    % searching deeps. 0 is default.
    [image_types_paths] = CheckImagePath(path_data,0);
    
    % all count info foe searched images at current folder.
    num_all_images = size(image_types_paths,1);
    num_L4_tm = 0;
    num_L5_tm = 0;
    num_L6_tm = 0;
    num_L7_tm_plus = 0;
    num_L8_oli_tirs = 0;
    num_S2A_msi = 0;
    num_S2B_msi = 0;
    
    % key info: the senor and path are needed as outputs, the first one is
    % to determine the default parameters for Fmask, and the second one is
    % to determine the path loading all data.
    sensors = [];
    paths = [];
    info_count_text = [];
    
    for i_all = 1:num_all_images
        % c: current.
        cimage_sensor = image_types_paths{i_all,1};
        cimage_num = image_types_paths{i_all,2};
        cimage_type = Convert2ImageType(cimage_sensor,cimage_num);
        clear cimage_num;
        cimage_path = image_types_paths{i_all,3};
        switch cimage_type
            case {'Landsat 4 TM'}
                num_L4_tm=num_L4_tm+1;
            case {'Landsat 5 TM'}
                num_L5_tm=num_L5_tm+1;
            case {'Landsat 6 TM'}
                num_L6_tm=num_L6_tm+1;
            case {'Landsat 7 ETM+'}
                num_L7_tm_plus=num_L7_tm_plus+1;
            case {'Landsat 8 OLI/TIRS'}
                num_L8_oli_tirs=num_L8_oli_tirs+1;
            case {'Sentinel 2A MSI'}
                % further check Sentinel 2 image folder, see there is
                % .SAFE.
                if contains(cimage_path,'.SAFE')
                    num_S2A_msi=num_S2A_msi+1;
                else
                    cimage_sensor = [];% no available Sentinel 2 data.
                end
            case {'Sentinel 2B MSI'}
                % see there is .SAFE.
                if contains(cimage_path,'.SAFE')
                    num_S2B_msi=num_S2B_msi+1;
                else
                    cimage_sensor = [];% no available Sentinel 2 data.
                end
        end
        if ~isempty(cimage_sensor)
            sensors = [sensors;{cimage_sensor}];
            paths = [paths;{cimage_path}];
        end
        clear cimage_sensor;
    end
    % renew num_all_images
    num_all_images = length(sensors);
    % used to notice user.
    text_line = 0;
    % multide
    if isequal(num_all_images,num_L4_tm)||...
            isequal(num_all_images,num_L5_tm)||...
            isequal(num_all_images,num_L6_tm)||...
            isequal(num_all_images,num_L7_tm_plus)||...
            isequal(num_all_images,num_L8_oli_tirs)||...
            isequal(num_all_images,num_S2A_msi)||...
            isequal(num_all_images,num_S2B_msi)
       % only for 1 type image
        if num_L4_tm > 0
        text_line = text_line+1;
        info_count_text{text_line} = sprintf('%s Landsat 4 TM images are found at ''%s''\n',...
            num2str(num_L4_tm),path_data);
        end
        if num_L5_tm > 0
            text_line = text_line+1;
            info_count_text{text_line} = sprintf('%s Landsat 5 TM images are found at ''%s''\n',...
                num2str(num_L5_tm),path_data);
        end
        if num_L6_tm > 0
            text_line = text_line+1;
            info_count_text{text_line} = sprintf('%s Landsat 6 TM images are found at ''%s''\n',...
            num2str(num_L6_tm),path_data);
        end
        if num_L7_tm_plus > 0
            text_line = text_line+1;
            info_count_text{text_line} = sprintf('%s Landsat 7 ETM+ images are found at ''%s''\n',...
                num2str(num_L7_tm_plus),path_data);
        end
        if num_L8_oli_tirs > 0
            text_line = text_line+1;
            info_count_text{text_line} = sprintf('%s Landsat 8 OLI/TIRS images are found at ''%s''\n',...
                 num2str(num_L8_oli_tirs),path_data);
        end
        if num_S2A_msi > 0
            text_line = text_line+1;
            info_count_text{text_line} = sprintf('%s Sentinel 2A MSI images are found at ''%s''\n',...
                num2str(num_S2A_msi),path_data);
        end
        if num_S2B_msi > 0
            text_line = text_line+1;
            info_count_text{text_line} = sprintf('%s Sentinel 2B MSI images are found at ''%s''\n',...
                num2str(num_S2B_msi),path_data);
        end
    else
         text_line = text_line+1;
         info_count_text{text_line} = sprintf(...
            'A total of %s images (as follows) are found at ''%s''\n',...
            num2str(num_all_images), path_data);
        if num_L4_tm > 0
        text_line = text_line+1;
        info_count_text{text_line} = sprintf('%s Landsat 4 TM images\n',...
            num2str(num_L4_tm));
        end
        if num_L5_tm > 0
            text_line = text_line+1;
            info_count_text{text_line} = sprintf('%s Landsat 5 TM images\n',...
                num2str(num_L5_tm));
        end
        if num_L6_tm > 0
            text_line = text_line+1;
            info_count_text{text_line} = sprintf('%s Landsat 6 TM images\n',...
            num2str(num_L6_tm));
        end
        if num_L7_tm_plus > 0
            text_line = text_line+1;
            info_count_text{text_line} = sprintf('%s Landsat 7 ETM+ images\n',...
                num2str(num_L7_tm_plus));
        end
        if num_L8_oli_tirs > 0
            text_line = text_line+1;
            info_count_text{text_line} = sprintf('%s Landsat 8 OLI/TIRS images\n',...
                 num2str(num_L8_oli_tirs));
        end
        if num_S2A_msi > 0
            text_line = text_line+1;
            info_count_text{text_line} = sprintf('%s Sentinel 2A MSI images\n',...
                num2str(num_S2A_msi));
        end
        if num_S2B_msi > 0
            text_line = text_line+1;
            info_count_text{text_line} = sprintf('%s Sentinel 2B MSI images\n',...
                num2str(num_S2B_msi));
        end
    end
    
    % if no data is found, give mention info.
    if isempty(sensors)
        text_line = text_line+1;
        info_count_text{text_line} = sprintf('%s available images are found at ''%s''\n',...
            '0',path_data);
        text_line = text_line+1;
        info_count_text{text_line} = sprintf('Please ensure the path is correct\n');
    end
%     fprintf(info_count);
end

function [image_types_paths] = CheckImagePath(path_data,subfolder_level)
%CHECKIMAGEPATH Ensure the input path is right to find a Landsats 4-8, and
%Sentinel 2 image(s).
% path_data - the input path
% subfolder_level - the level of subfolders that can be used to limit
% searching deeps. 0 is default.
    
    % If the searching deeps are more than 5, stop and return;
    if subfolder_level > 5
        image_types_paths = [];
        return;
    end
	image_types_paths = [];
    % first search the image(s) at current folder.image_types_paths
    [sensor,num_Lst,~,~] = LoadSensorType(path_data);
    if isempty(sensor)
        % if no available image at current folder,
        % and search its subfolders.
          subfolders = dir(path_data);
          for i_sub=1:length(subfolders)
              % filter out the file names starting with '.', that is not
              % right folder (system crashes).
              if strcmp(subfolders(i_sub).name(1),'.')
                  continue;
              end
              % go to search the images at each subfolder
              path_subfoler = fullfile(subfolders(i_sub).folder,...
                  subfolders(i_sub).name);
              [image_types_paths_sub] = CheckImagePath(path_subfoler,subfolder_level+1);
              if ~isempty(image_types_paths_sub)
                  image_types_paths =[image_types_paths;image_types_paths_sub];
              end
          end
          
    else
        % successfully searched a supported image.
        % and, return the sensor and the image path
        image_types_paths = [image_types_paths;{sensor,num_Lst,path_data}];
    end
end
function image_type = Convert2ImageType(sensor,num_Lst)
% CONVERT@IMAGETYPE Contruct image type from the sensor and num

    switch sensor
        case 'L_TM'
            image_type = ['Landsat ',num_Lst,' TM'];
        case 'L_ETM_PLUS'
            image_type = ['Landsat ',num_Lst,' ETM+'];
        case 'L_OLI_TIRS'
            image_type = ['Landsat ',num_Lst,' OLI/TIRS'];
        case 'S_MSI'
            image_type = ['Sentinel ',num_Lst,' MSI'];
        otherwise
            image_type=[];
%             error(['Errors occur when searching images at ''', path_data],'''.');
    end
end

