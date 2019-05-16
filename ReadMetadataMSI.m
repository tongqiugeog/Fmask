function MSIinfo=ReadMetadataMSI(FilesInfo)
% DataStrip='S2A_OPER_PRD_MSIL1C_PDMC_20150818T101237_R022_V20150813T102406_20150813T102406';
% Granule='S2A_OPER_MSI_L1C_TL_MTI__20150813T201603_A000734_T32TPS_N01';

DirIn=FilesInfo.DirIn;
% Num=FilesInfo.Num;
DataStrip = FilesInfo.DataStrip; 
Granule=FilesInfo.Granule;

Dmain = DirIn ;

if Granule(1) == 'L'
    FileName = fullfile(Dmain, [DataStrip '.SAFE'],'GRANULE', Granule, 'MTD_TL.xml');
%     FileName = [Dmain DataStrip '.SAFE/GRANULE/' Granule '/' 'MTD_TL.xml'];
else
%     FileName = [Dmain DataStrip '.SAFE/GRANULE/' Granule '/' strrep(Granule(1:55),'MSI','MTD') '.xml'];
    FileName = fullfile(Dmain,[DataStrip '.SAFE'],'GRANULE', Granule,[strrep(Granule(1:55),'MSI','MTD') '.xml']);
end
    
S = xml2struct(FileName);

S = S.n1_colon_Level_dash_1C_Tile_ID;

TileID = S.n1_colon_General_Info.TILE_ID.Text;

ds =  S.n1_colon_General_Info.SENSING_TIME.Text ;
DN = (datenum( ds ,'yyyy-mm-ddTHH:MM:SS')+str2double(ds(end-4:end-1))/3600/24);

DS = datestr(DN);

bandIdList = {'01' '02' '03' '04' '05' '06' '07' '08' '8A' '09' '10' '11' '12'} ;

GeoInfo.UTM_zone = S.n1_colon_Geometric_Info.Tile_Geocoding.HORIZONTAL_CS_NAME.Text(end-2:end);

%%
M = S.n1_colon_Quality_Indicators_Info.Pixel_Level_QI.MASK_FILENAME ;

for i=1:length(M)
  if isfield(M{i}.Attributes,'bandId')
  Mask.(M{i}.Attributes.type).(['B' M{i}.Attributes.bandId]) = M{i}.Text ;
  else
    Mask.(M{i}.Attributes.type) = M{i}.Text ;
  end
end
%%
GeoInfo.Size.R10 = [str2num(S.n1_colon_Geometric_Info.Tile_Geocoding.Size{1}.NCOLS.Text) str2num(S.n1_colon_Geometric_Info.Tile_Geocoding.Size{1}.NROWS.Text)] ;
GeoInfo.Size.R20 = [str2num(S.n1_colon_Geometric_Info.Tile_Geocoding.Size{2}.NCOLS.Text) str2num(S.n1_colon_Geometric_Info.Tile_Geocoding.Size{2}.NROWS.Text)] ;
GeoInfo.Size.R60 = [str2num(S.n1_colon_Geometric_Info.Tile_Geocoding.Size{3}.NCOLS.Text) str2num(S.n1_colon_Geometric_Info.Tile_Geocoding.Size{3}.NROWS.Text)] ;

GeoInfo.Xstart.R10 =  str2num(S.n1_colon_Geometric_Info.Tile_Geocoding.Geoposition{1}.ULX.Text);
GeoInfo.Xstart.R20 =  str2num(S.n1_colon_Geometric_Info.Tile_Geocoding.Geoposition{2}.ULX.Text);
GeoInfo.Xstart.R60 =  str2num(S.n1_colon_Geometric_Info.Tile_Geocoding.Geoposition{3}.ULX.Text);

GeoInfo.Ystart.R10 =  str2num(S.n1_colon_Geometric_Info.Tile_Geocoding.Geoposition{1}.ULY.Text);
GeoInfo.Ystart.R20 =  str2num(S.n1_colon_Geometric_Info.Tile_Geocoding.Geoposition{2}.ULY.Text);
GeoInfo.Ystart.R60 =  str2num(S.n1_colon_Geometric_Info.Tile_Geocoding.Geoposition{3}.ULY.Text);

%%
clear dum
for i=1:23
  dum(i,:) = str2num(S.n1_colon_Geometric_Info.Tile_Angles.Sun_Angles_Grid.Zenith.Values_List.VALUES{i}.Text) ;
end
Angles.SZA = dum ;

clear dum
for i=1:length(S.n1_colon_Geometric_Info.Tile_Angles.Sun_Angles_Grid.Azimuth.Values_List.VALUES)
  dum(i,:) = str2num(S.n1_colon_Geometric_Info.Tile_Angles.Sun_Angles_Grid.Azimuth.Values_List.VALUES{i}.Text) ;
end
Angles.SAA = dum ;

bandIdList = 0:12 ;
DeteIdList = 1:12 ;

Angles.VZA = nan(23,23,length(bandIdList),length(DeteIdList));
Angles.VAA = nan(23,23,length(bandIdList),length(DeteIdList));

Sdum = S.n1_colon_Geometric_Info.Tile_Angles.Viewing_Incidence_Angles_Grids;
for j=1:length(Sdum)
  bandId = str2num(Sdum{j}.Attributes.bandId) ;
  DeteId = str2num(Sdum{j}.Attributes.detectorId) ;
  clear dum
  for i=1:23
    dum(i,:) = str2num(Sdum{j}.Zenith.Values_List.VALUES{i}.Text) ;
  end
  Angles.VZA(:,:,bandIdList==bandId,DeteIdList==DeteId) = dum ;
  clear dum
  for i=1:23
    dum(i,:) = str2num(Sdum{j}.Azimuth.Values_List.VALUES{i}.Text) ;
  end
  Angles.VAA(:,:,bandIdList==bandId,DeteIdList==DeteId) = dum ;
end

Angles.Mean.SZA = str2num(S.n1_colon_Geometric_Info.Tile_Angles.Mean_Sun_Angle.ZENITH_ANGLE.Text) ;
Angles.Mean.SAA = str2num(S.n1_colon_Geometric_Info.Tile_Angles.Mean_Sun_Angle.AZIMUTH_ANGLE.Text) ;

Angles.Mean.VZA2 = nan(length(bandIdList),1);
Angles.Mean.VAA2 = nan(length(bandIdList),1);

Sdum = S.n1_colon_Geometric_Info.Tile_Angles.Mean_Viewing_Incidence_Angle_List.Mean_Viewing_Incidence_Angle;
for j=1:length(Sdum)
  bandId = str2num(Sdum{j}.Attributes.bandId) ;
  Angles.Mean.VZA2(bandIdList==bandId) = str2num(Sdum{j}.ZENITH_ANGLE.Text) ;
  Angles.Mean.VAA2(bandIdList==bandId) = str2num(Sdum{j}.AZIMUTH_ANGLE.Text) ;
end

Angles.Mean.VZA = mean(Angles.Mean.VZA2);
Angles.Mean.VAA = mean(Angles.Mean.VAA2);

Angles.bandIdList=bandIdList;
Angles.DeteIdList=DeteIdList;
%%
MSIinfo.TileID = TileID;
MSIinfo.DN = DN;
MSIinfo.DS = DS;
MSIinfo.bandIdList = bandIdList ;
MSIinfo.Angles = Angles;
MSIinfo.GeoInfo = GeoInfo;
MSIinfo.Mask = Mask ;
