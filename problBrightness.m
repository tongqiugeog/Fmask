function prob_brightness = problBrightness(HOT,idlnd,l_pt,h_pt)
%PROBLBRIGHTNESS calculate brightness probability using HOT

    F_hot=HOT(idlnd); % get clear HOTs
    clear idlnd;
    % 0.175 percentile background HOT (low)
    t_hotl=prctile(F_hot,100*l_pt)-400;
    clear l_pt;
    % 0.825 percentile background HOT (high)
    t_hoth=prctile(F_hot,100*h_pt)+400;
    clear F_hot h_pt;

    prob_brightness=(HOT-t_hotl)./(t_hoth-t_hotl);
    clear HOT t_hotl t_hoth;
    prob_brightness(prob_brightness<0)=0;
    prob_brightness(prob_brightness>1)=1; % this cannot be higher 1 (maybe have commission errors from bright surfaces).
end

