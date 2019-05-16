function [MSIinfo Angles CM]  = ReadSunViewGeometryMSI (DataStrip,Granule,BandSel,PsizeOut,Dmain)
% DataStrip='S2A_OPER_PRD_MSIL1C_PDMC_20150818T101451_R080_V20150817T114755_20150817T114755';
% Granule='S2A_OPER_MSI_L1C_TL_SGS__20150817T131818_A000792_T28QBG_N01';
% this is derived from Fmask 3.3 for Sentinel 2.

InterpMethod = 'linear';'nearest';

% Dmain = '/nobackupp6/jju/S2SAFE/';
Pppplot = 0 ;
% PsizeOut = 20 ;
if nargin<5
  FilesInfo.DirIn = CodeStartingLinesSentinel(-1);
else
  FilesInfo.DirIn = Dmain;
end
% FilesInfo.Num = '.04';
FilesInfo.DataStrip = DataStrip;
FilesInfo.Granule=Granule;

if Granule(1) == 'L'
    mainpath = fullfile(Dmain,[DataStrip '.SAFE']);
else
    mainpath = fullfile(Dmain, [DataStrip '.SAFE'],'GRANULE', Granule, 'QI_DATA');
end

try
    MSIinfo=ReadMetadataMSI(FilesInfo);
catch
    error('The metadata was wrongly loaded. Please input Sentinel-2 TOA data.');
end

if nargout==1
  return
end
bandIdList = {'01' '02' '03' '04' '05' '06' '07' '08' '8A' '09' '10' '11' '12'} ;
BandList = 0:12 ;
if exist('BandSel','var')
  bandIdList = bandIdList(BandSel);
  BandList = BandList(BandSel);
end
if ~exist('PsizeOut','var')
  PsizeOut = 10 ;
end
% bandIdList = bandIdList (8:9);
% BandList = BandList (8:9);


% MethInterp100m10m = 'nearest';

% [BandList]=find(MSIinfo.Angles.Mean.VZA2==median(MSIinfo.Angles.Mean.VZA2)) ;

%%
X5km = MSIinfo.GeoInfo.Xstart.R10 : 5000 : MSIinfo.GeoInfo.Xstart.R10 + 5000*22 ;
Y5km = MSIinfo.GeoInfo.Ystart.R10 :-5000 : MSIinfo.GeoInfo.Ystart.R10 - 5000*22 ;
% [X5kmat Y5kmat]=meshgrid(X5km,Y5km);
% Psize = 100 ;
% X100m = single(MSIinfo.GeoInfo.Xstart.R20+Psize/2 : Psize : MSIinfo.GeoInfo.Xstart.R20-Psize/2 + Psize*(MSIinfo.GeoInfo.Size.R10(2)/10))';
% Y100m = single(MSIinfo.GeoInfo.Ystart.R20-Psize/2 :-Psize : MSIinfo.GeoInfo.Ystart.R20+Psize/2 - Psize*(MSIinfo.GeoInfo.Size.R10(1)/10))';
% [X100mat Y100mat]=meshgrid(X100m,Y100m);

Psize = 10 ;
X10m = single(MSIinfo.GeoInfo.Xstart.R20 : Psize : MSIinfo.GeoInfo.Xstart.R20+(MSIinfo.GeoInfo.Size.R10(2)-1)*Psize)  + Psize/2  ;
Y10m = single(MSIinfo.GeoInfo.Ystart.R20 :-Psize : MSIinfo.GeoInfo.Ystart.R20-(MSIinfo.GeoInfo.Size.R10(1)-1)*Psize)' - Psize/2  ;

% X100m = imresize(X10m,10/100,'box') ;
% Y100m = imresize(Y10m,10/100,'box') ;

XPsizeOut = imresize(X10m,10/PsizeOut,'box') ;
YPsizeOut = imresize(Y10m,10/PsizeOut,'box') ;
clear PsizeOut;

% Psize = PsizeOut ;

% [X100mat Y100mat]=meshgrid(X100m,Y100m);
[XPsizeOutmat YPsizeOutmat]=meshgrid(XPsizeOut,YPsizeOut);


if Pppplot == 1
  %   Psize = PsizeOut ;
  %   X10m = single(MSIinfo.GeoInfo.Xstart.R20+Psize/2 : Psize : MSIinfo.GeoInfo.Xstart.R20-Psize/2 + (MSIinfo.GeoInfo.Size.R10(2)*10/Psize))';
  %   Y10m = single(MSIinfo.GeoInfo.Ystart.R20-Psize/2 :-Psize : MSIinfo.GeoInfo.Ystart.R20+Psize/2 - (MSIinfo.GeoInfo.Size.R10(1)*10/Psize))';
  %   [X10mat Y100mat]=meshgrid(X100m,Y100m);
  [X5kmat Y5kmat]=meshgrid(X5km,Y5km);
end
% % kern = ones(20);imcircle(20);

%% cloud Mask

if nargout==3
  
  FileName = fullfile(mainpath, MSIinfo.Mask.MSK_CLOUDS);
  S = xml2struct(FileName);
  
  if isfield(S.eop_colon_Mask,'eop_colon_maskMembers')
    
    S = S.eop_colon_Mask.eop_colon_maskMembers;
    
    
    
    
    
    if length(S)==1
      %     DetectFootPrint.Ndect(1) = str2num(S(1).eop_colon_MaskFeature.Attributes.gml_colon_id(24:25));
      %     disp('attention')
      %     DetectFootPrint.Nband(1) = find(strcmp(MSIinfo.bandIdList , S(1).eop_colon_MaskFeature.Attributes.gml_colon_id(21:22)));
      Vect = str2num(S(1).eop_colon_MaskFeature.eop_colon_extentOf.gml_colon_Polygon.gml_colon_exterior.gml_colon_LinearRing.gml_colon_posList.Text);
      Vect = reshape(Vect' , [2 length(Vect)/2 ])';
      CloudMask.Data{1,1} = Vect ;
      %     Z = MSIinfo.Angles.VZA(:,:,MSIinfo.Angles.bandIdList==DetectFootPrint.Nband(1),MSIinfo.Angles.DeteIdList==DetectFootPrint.Ndect(1)) ;
    else
      
      for i = 1:length(S)
        %       DetectFootPrint.Ndect(i,1) = str2num(S{i}.eop_colon_MaskFeature.Attributes.gml_colon_id(24:25));
        %       DetectFootPrint.Nband(i,1) = strcmp(MSIinfo.bandIdList , S{i}.eop_colon_MaskFeature.Attributes.gml_colon_id(21:22));
        Vect = str2num(S{i}.eop_colon_MaskFeature.eop_colon_extentOf.gml_colon_Polygon.gml_colon_exterior.gml_colon_LinearRing.gml_colon_posList.Text);
        Vect = reshape(Vect' , [2 length(Vect)/2 ])';
        CloudMask.Data{i,1} = Vect ;
        %       if Pppplot == 1
        %         IDX = knnsearch(X10m,DetectFootPrint.Data{i,1}(:,1));
        %         Ximg = X10m(IDX);
        %         IDX = knnsearch(Y10m,DetectFootPrint.Data{i,1}(:,2));
        %         Yimg = Y10m(IDX);
        %         plot(Ximg,Yimg,'k')
        %         Z = MSIinfo.Angles.VAA(:,:,MSIinfo.Angles.bandIdList==DetectFootPrint.Nband(i),MSIinfo.Angles.DeteIdList==DetectFootPrint.Ndect(i)) ;
        %         scatter(X5kmat(:),Y5kmat(:),30,Myreshape(Z),sym(i))
        %       end
        
      end
    end
    
    Vect = [NaN NaN];
    for i = 1:length(S)
      Vect = cat(1,Vect ,CloudMask.Data{i},[NaN NaN]);
    end
    
    %   X10m = single(MSIinfo.GeoInfo.Xstart.R20+Psize/2 : Psize : MSIinfo.GeoInfo.Xstart.R10-Psize/2 + Psize*(MSIinfo.GeoInfo.Size.R10(2)))';
    %   Y10m = single(MSIinfo.GeoInfo.Ystart.R20-Psize/2 :-Psize : MSIinfo.GeoInfo.Ystart.R10+Psize/2 - Psize*(MSIinfo.GeoInfo.Size.R10(1)))';
    
    IDX = knnsearch(XPsizeOut',Vect(:,1));
    IDY = knnsearch(YPsizeOut,Vect(:,2));
    CM = uint8(poly2mask(double(IDX),double(IDY),length(YPsizeOut),length(XPsizeOut)));
	clear IDX IDY;
  else
    CM = zeros( length(YPsizeOut),length(XPsizeOut) ,'uint8');
  end
end


%%
%%
%%
%%
%% Angles
for iB = 1:length(BandList)
  %   tic
  FileName = fullfile(mainpath,...
    MSIinfo.Mask.MSK_DETFOO.(['B' num2str(BandList(iB))]));
  
  S = xml2struct(FileName);
  clear FileName;
  
  if isfield(S.eop_colon_Mask.eop_colon_maskMembers,'eop_colon_MaskFeature')
    S = S.eop_colon_Mask.eop_colon_maskMembers.eop_colon_MaskFeature;
  else
    S = S.eop_colon_Mask.eop_colon_maskMembers;
  end
  %   figure
  
  if length(S)==1

    try
        % the following codes are unworkable anymore.
        if ~(isfield(S,'eop_colon_MaskFeature'))
             S{1}.eop_colon_MaskFeature = S{1};
        end
        DetectFootPrint.Ndect(1) = str2num(S{1}.eop_colon_MaskFeature.Attributes.gml_colon_id(24:25));
        disp('attention')
        DetectFootPrint.Nband(1) = find(strcmp(MSIinfo.bandIdList , S(1).eop_colon_MaskFeature.Attributes.gml_colon_id(21:22)));
        Vect = str2num(S(1).eop_colon_MaskFeature.eop_colon_extentOf.gml_colon_Polygon.gml_colon_exterior.gml_colon_LinearRing.gml_colon_posList.Text);
        Vect = reshape(Vect' , [3 length(Vect)/3 ])';
        DetectFootPrint.Data{1,1} = Vect ;
        clear Vect;
    catch
        if ~(isfield(S,'eop_colon_MaskFeature'))
            S.eop_colon_MaskFeature = S;
        end
        DetectFootPrint.Ndect(1) = str2num(S.eop_colon_MaskFeature.Attributes.gml_colon_id(24:25));
        disp('attention')
        DetectFootPrint.Nband(1) = find(strcmp(bandIdList , S.eop_colon_MaskFeature.Attributes.gml_colon_id(21:22)));
        Vect = str2num(S.eop_colon_MaskFeature.eop_colon_extentOf.gml_colon_Polygon.gml_colon_exterior.gml_colon_LinearRing.gml_colon_posList.Text);
        Vect = reshape(Vect' , [3 length(Vect)/3 ])';
        DetectFootPrint.Data{1,1} = Vect ;
        clear Vect;
    end
    %     Z = MSIinfo.Angles.VZA(:,:,MSIinfo.Angles.bandIdList==DetectFootPrint.Nband(1),MSIinfo.Angles.DeteIdList==DetectFootPrint.Ndect(1)) ;
  else
    if Pppplot == 1
      figure
      hold on
      sym = 'oxoxoxoxox';
    end
    for i = 1:length(S)
      if ~(isfield(S{i},'eop_colon_MaskFeature'))
        S{i}.eop_colon_MaskFeature = S{i};
      end
      DetectFootPrint.Ndect(i,1) = str2num(S{i}.eop_colon_MaskFeature.Attributes.gml_colon_id(24:25));
      %       DetectFootPrint.Nband(i,1) = strcmp(MSIinfo.bandIdList , S{i}.eop_colon_MaskFeature.Attributes.gml_colon_id(21:22));
      DetectFootPrint.Nband(i,1) = strcmp(bandIdList , S{i}.eop_colon_MaskFeature.Attributes.gml_colon_id(21:22));
      
      Vect = str2num(S{i}.eop_colon_MaskFeature.eop_colon_extentOf.gml_colon_Polygon.gml_colon_exterior.gml_colon_LinearRing.gml_colon_posList.Text);
      Vect = reshape(Vect' , [3 length(Vect)/3 ])';
      DetectFootPrint.Data{i,1} = Vect ;
      clear Vect;
      if Pppplot == 1
        
        IDX = knnsearch(X10m,DetectFootPrint.Data{i,1}(:,1));
        Ximg = X10m(IDX);
        IDX = knnsearch(Y10m,DetectFootPrint.Data{i,1}(:,2));
        Yimg = Y10m(IDX);
        clear IDX;
        plot(Ximg,Yimg,'k')
        clear Ximg Yimg;
        Z = MSIinfo.Angles.VAA(:,:,MSIinfo.Angles.bandIdList==DetectFootPrint.Nband(i),MSIinfo.Angles.DeteIdList==DetectFootPrint.Ndect(i)) ;
        scatter(X5kmat(:),Y5kmat(:),30,Myreshape(Z),sym(i))
      end
    end
  end
  clear S;
  %   %   axis image
  
  
  %% 5km => 10m
  Matdetec = findMatdetecFootprint(DetectFootPrint,XPsizeOut,YPsizeOut);
  
  %   JJ.B04_detfoo = hdfread('D:\UMD\Data\HLS\Products\detfoo.T36JTT.2015350.v1.0.hdf', '/B04_detfoo', 'Index', {[1  1],[1  1],[1830  1830]});
  
  Angles.Matdetec(:,:,iB) = uint8(Matdetec) ;
  Angles.DetectFootPrint{iB} = DetectFootPrint ;
  % return
  clear dum*
  
  % figure
  % hold on
  % imagesc(X10m,Y10m,Matdetec)
  % for i = 1:length(S)
  %   plot(DetectFootPrint.Data{i,1}(:,1),DetectFootPrint.Data{i,1}(:,2),'k')
  % end
  
  Chp = {'VAA' 'VZA'} ;
  
  for iC=1:2
    Angles.([Chp{iC} '_B' bandIdList{iB}]) = nan(size(Matdetec),'single');
    for i = 1:length(DetectFootPrint.Ndect)
      
      Z = MSIinfo.Angles.(Chp{iC})(:,:,BandSel,MSIinfo.Angles.DeteIdList==DetectFootPrint.Ndect(i)) ;
      Z = single(inpaint_nans(Z));
      
      testMatdetec = Matdetec==i ;
      [Ista dum] = find(max(testMatdetec,[],2),1,'first');
      [Iend dum] = find(max(testMatdetec,[],2),1,'last');
      [dum Jsta] = find(max(testMatdetec,[],1),1,'first');
      [dum Jend] = find(max(testMatdetec,[],1),1,'last');
      clear dum;
      %       Ista = floor(Ista/10)+1;
      %       Jsta = floor(Jsta/10)+1;
      %       Iend = ceil(Iend/10);
      %       Jend = ceil(Jend/10);
      
      ZI = zeros(size(Matdetec),'single');
      
      ZI(Ista:Iend,Jsta:Jend) = interp2(X5km,Y5km,Z,XPsizeOutmat(Ista:Iend,Jsta:Jend),YPsizeOutmat(Ista:Iend,Jsta:Jend),InterpMethod);
      clear Z Ista Iend Jsta Jend;
      %       ZI(Ista:Iend,Jsta:Jend) = interp2(X5km,Y5km,Z,X10mat(Ista:Iend,Jsta:Jend),Y100mat(Ista:Iend,Jsta:Jend),InterpMethod);
      %       ZI = imresize(ZI,10,MethInterp100m10m);
      
      Angles.([Chp{iC} '_B' bandIdList{iB}])(testMatdetec) = ZI(testMatdetec) ;
      clear ZI testMatdetec;
    end
  end
  clear Matdetec;
  %   toc
end


Z = MSIinfo.Angles.SZA ;
ZI = interp2(X5km,Y5km,Z,XPsizeOutmat,YPsizeOutmat,InterpMethod);
clear Z;
% ZI=imresize(ZI,10,MethInterp100m10m);
Angles.SZA = ZI;

Z = MSIinfo.Angles.SAA ;
ZI = interp2(X5km,Y5km,Z,XPsizeOutmat,YPsizeOutmat,InterpMethod);
clear Z X5km Y5km XPsizeOutmat YPsizeOutmat InterpMethod;
% ZI=imresize(ZI,10,MethInterp100m10m);
Angles.SAA = ZI;
clear ZId;

% for iB = BandList
%   Angles.(['RAA_B' bandIdList{iB+1}]) = abs(Angles.(['VAA_B' num2str(iB)]) - Angles.SAA) ;
%   Angles=rmfield(Angles,['VAA_B' num2str(iB)]);
% end
% Angles=rmfield(Angles,'SAA');

fn = fieldnames(Angles);
dum = cellfun(@(x) isempty(strfind(x,'A')),fn);
fn(dum) = [] ;
clear dum;
for i=1:length(fn)-1
  dum = Angles.(fn{i}) ;
  dum(dum<0) = 655.36 ;
  Angles.(fn{i}) = uint16(dum*100);
  clear dum;
end
clear fn;

Angles.BandList = BandList ;
Angles.bandIdList = bandIdList ;
clear BandList bandIdList;

return

y = structfun(@(x) max(x(:)), Angles)
