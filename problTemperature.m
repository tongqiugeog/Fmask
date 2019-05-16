function [prob_temp, t_templ,t_temph] = problTemperature(BandBT,idclr,l_pt,h_pt)
%PROBTEMPERATURE calculate temperature probability for land respectively.

    % [Temperature test (over land)]
    F_temp=BandBT(idclr); % get clear temperature
    clear idclr;
    t_buffer=4*100;
    % 0.175 percentile background temperature (low)
    t_templ=prctile(F_temp,100*l_pt);
    % 0.825 percentile background temperature (high)
    t_temph=prctile(F_temp,100*h_pt);
    clear F_temp l_pt h_pt;

    t_tempL=t_templ-t_buffer;
    t_tempH=t_temph+t_buffer;
    clear t_buffer;
    Temp_l=t_tempH-t_tempL;
    clear t_tempL;
    prob_temp=(t_tempH-BandBT)/Temp_l;
    clear BandBT t_tempH Temp_l;
    % Temperature can have prob > 1
    prob_temp(prob_temp<0)=0;
%     prob_temp(prob_temp>1)=1;
    
end

