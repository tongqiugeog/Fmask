function prob_vari = problSpectralVaribility(ndvi,ndsi,ndbi,whiteness,SatuGreen,SatuRed)
%PROBLSPECTRALVARIBILITY 
    % [varibility test over land]
    ndsi(SatuGreen==true&ndsi<0)=0;
    clear SatuGreen;
    ndvi(SatuRed==true&ndvi>0)=0;
    clear SatuRed;
    prob_vari=1-max(max(max(abs(ndsi),abs(ndvi)),abs(ndbi)),whiteness);
    clear ndsi ndvi ndbi whiteness;
end

