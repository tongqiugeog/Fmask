classdef FmaskParameters
    %FMASKPARAMETERS Save all fmask parameters here.
    
    properties
        CloudBuffer % pixels
        CloudShadowBuffer % pixels
        SnowBuffer % pixels
        ThinWeight
        CloudProbabilityThershold
        PFPCErosionRadius  % meters
        OutputResolution  % meters
        ShadowWater % yes or no for masking cloud shadow over water
    end
    
    methods
        function obj = FmaskParameters(sensor)
            %FMASKPARAMETERS Construct an instance of this class according
            %to input image.
            
            % public and constant paras.
            obj.CloudBuffer=3;
            obj.CloudShadowBuffer=3;
            obj.SnowBuffer=0;
            
            % mask out the shadow of the cloud over water?
            % default: we do not provide the cloud shadow over water since
            % this processing will be very time-comsuing but less meanful.
            obj.ShadowWater=0; 
            
            % different paras for different sensors.
            switch sensor
                case 'S_MSI'
                    obj.ThinWeight=0.5;
                    obj.CloudProbabilityThershold=20.00;
                    obj.OutputResolution=20;
                    obj.PFPCErosionRadius=90;% mirrored from Landsat 8.
                case 'L_OLI_TIRS'
                    obj.ThinWeight=0.3;
                    obj.CloudProbabilityThershold=17.50;
                    obj.OutputResolution=30;
                    obj.PFPCErosionRadius=90;
                case {'L_TM','L_ETM_PLUS'}
                    obj.ThinWeight=0.0;
                    obj.CloudProbabilityThershold=10.00;
                    obj.OutputResolution=30;
                    obj.PFPCErosionRadius=150;
            end
        end
    end
end

