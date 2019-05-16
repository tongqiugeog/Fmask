classdef ObjTOABT
    %OBJTOA TOA & BT data
%     Top of Atmosphere reflectance (TOA) and Brightness Temperature (BT)
%     for all used bands. * indicates the used bands.
%     ------------------------------------------------------------------------------------------------------------
%                Satellite                             Landsat                                     Sentinel
%     Band Name  -------------------------------------------------------------------------------------------------
%                Sensor                TM                ETM+                OLI/TIRS                 MSI
%                -------------------------------------------------------------------------------------------------
%     CA                              NULL               NULL             1:0.43-0.45(30m)      1:0.433-0.453(60m)
%     Blue*                     1:0.45-0.515(30m)   1:0.45-0.515(30m)     2:0.452-0.512(30m)    2:0.458-0.523(10m)
%     Green*                    2:0.52-0.60(30m)    2:0.52-0.60(30m)      3:0.533-0.590(30m)    3:0.543-0.578(10m)
%     Red*                      3:0.63-0.69(30m)    3:0.63-0.69(30m)      4:0.636-0.673(30m)    4:0.650-0.680(10m)
%     VRE1                            NULL               NULL                   NULL            5:0.698-0.713(20m)
%     VRE2                            NULL               NULL                   NULL            6:0.733-0.748(20m)
%     VRE3                            NULL               NULL                   NULL            7:0.765-0.785(20m)
%     NIR*                      4:0.75-0.90(30m)    4:0.75-0.90(30m)      5:0.851-0.879(30m)    8:0.855-0.875(10m)
%     NNIR                            NULL               NULL                   NULL            8A(20m)
%     WV                              NULL               NULL                   NULL            9:0.930-0.950(60m) 
%     Cirrus*                         NULL               NULL             9:1.360–1.390(30m)   10:1.365-1.385(60m)
%     SWIR1*                    5:1.55-1.75(30m)    5:1.55-1.75(30m)      6:1.566-1.651(30m)   11:1.565-1.655(20m)
%     SWIR2*                    7:2.09-2.35(30m)    7:2.09-2.35(30m)      7:2.107-2.294(30m)   12:2.100-2.280(20m)
%     Pan                       8:0.52-0.90(15m)    8:0.500-0.680(15m)          NULL                  NULL
%     BT                        6:10.40-12.50(60m)  6:10.40-12.50(60m)   10:10.6-11.2(100m)           NULL
%     ------------------------------------------------------------------------------------------------------------
    properties
        %% TOA
        BandCA % Coastal Aerosol (for S2)
        BandBlue
        BandGreen
        BandRed
        BandVRE1 %Vegetation Red Edge at 0.705 (for S2)
        BandVRE2 %Vegetation Red Edge at 0.740 (for S2)
        BandVRE3 %Vegetation Red Edge at 0.783 (for S2) band 7
        BandNIR8 % Narrow NIR (for S2) band 8
        BandNIR % Narrow NNIR (for S2) band 8a
        BandWV % Water Vapour
        BandCirrus % (for L8 and S2)
        BandSWIR1
        BandSWIR2
        %% BT
        BandBT %Brightness Temperature (Unit: Celsius degree) (for L4-8) 
        %% Statured Visible Bands
        SatuBlue
        SatuGreen
        SatuRed
    end
    
    methods
    end
    
end

