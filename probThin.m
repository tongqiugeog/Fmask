function prob_thin = probThin(BandCirrus)
%PROBCIRRUS
    prob_thin = BandCirrus./400;
    clear BandCirrus;
%     prob_thin(prob_thin>1)=1;
    prob_thin(prob_thin<0)=0;
end

