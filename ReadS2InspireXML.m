function bbox = ReadS2InspireXML(xml_path)
%READS2INSPIREXML Load INSPIRE XML data at .SAFE directory, and can give us
%image extent.
% The range of Sentinel2's longtitude is [0 360] sometimes, which need to be [-180 180]  by Shi 3/17/2020
    FileName = fullfile(xml_path,'INSPIRE.xml');
    S = xml2struct(FileName);
    clear FileName;
    S = S.gmd_colon_MD_Metadata.gmd_colon_identificationInfo.gmd_colon_MD_DataIdentification;
    extentS=S.gmd_colon_extent.gmd_colon_EX_Extent.gmd_colon_geographicElement.gmd_colon_EX_GeographicBoundingBox;
    clear S;
    % bbox=[north,south,west,east];
    bbox=[str2double(extentS.gmd_colon_northBoundLatitude.gco_colon_Decimal.Text),...
        str2double(extentS.gmd_colon_southBoundLatitude.gco_colon_Decimal.Text),...
        str2double(extentS.gmd_colon_westBoundLongitude.gco_colon_Decimal.Text),...
        str2double(extentS.gmd_colon_eastBoundLongitude.gco_colon_Decimal.Text)];
    
    if bbox(3) > 180
        bbox(3) = bbox(3) - 360;
    end
    if bbox(4) > 180
        bbox(4) = bbox(4) - 360;
    end
end

