function wprob_brightness = probwBrightness(BandSWIR)
%PROBWBRIGHTNESS calculate brightness probability
 % [Brightness test (over water)]
    t_bright=1100;
    wprob_brightness=BandSWIR./t_bright;
    clear BandSWIR t_bright;
    wprob_brightness(wprob_brightness>1)=1;
    wprob_brightness(wprob_brightness<0)=0;
end

