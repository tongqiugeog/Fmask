function satu_Bv = Saturate( satu_blue,satu_green,satu_red )
%SATURE Summary of this function goes here
%   Detailed explanation goes here
    satu_Bv=satu_blue+satu_green+satu_red>=1; % saturated pixes at any visible bands.
    clear satu_blue satu_green satu_red;
end

