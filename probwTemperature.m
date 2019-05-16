function wprob_temp = probwTemperature(BandBT,idwt,h_pt)
%PROBTEMPERATURE calculate temperature probability for water.
    
    % [temperature test (over water)]
    F_wtemp=BandBT(idwt); % get clear water temperature
    clear idwt;
    t_wtemp=prctile(F_wtemp,100*h_pt);
    clear F_wtemp h_pt;
    wprob_temp=(t_wtemp-BandBT)/400;
    clear t_wtemp BandBT;
    wprob_temp(wprob_temp<0)=0;
    %    wTemp_prob(wTemp_prob > 1) = 1;
end

