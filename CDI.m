function cdi = CDI(S2band7,S2band8,S2band8A)
%CDI This is used seprate bright surface from cloud by following David,
%2018 RSE
    if ~isempty(S2band7)&&~isempty(S2band8)
        ratio_8A_8 = S2band8./S2band8A;
        clear S2band8;
        ratio_8A_7 = S2band7./S2band8A;
        clear S2band7 S2band8A;
        
        std_ratio_8A_8 =stdfilt(ratio_8A_8, true(7));
        clear ratio_8A_8;
        var_ratio_8A_8 = std_ratio_8A_8 .^2;
        clear std_ratio_8A_8;
        
        std_ratio_8A_7 = stdfilt(ratio_8A_7, true(7));
        clear ratio_8A_7;
        var_ratio_8A_7 = std_ratio_8A_7 .^2;
        clear std_ratio_8A_7;

        cdi = (var_ratio_8A_7-var_ratio_8A_8)./(var_ratio_8A_8+var_ratio_8A_7);
        
%         dim = size(cdi);
%         cdi = imresize(cdi,1/7,'box');
%         cdi = ordfilt2(cdi,1,ones(3,3));
%         cdi = imresize(cdi,dim,'nearest');
%         clear dim;
        
        
%         h = fspecial('average', 21);
%         cdi = filter2(h, cdi);
    else
        cdi = [];
    end
end


