function line_enhanced = EnhanceLine(band)
%ENHANCELINE enhance NDBI to light urban/built-up areas and dark other
%bright surface, such as desert, rock. ref Guindon et. RSE 2004

% also can see the details at
% https://homepages.inf.ed.ac.uk/rbf/HIPR2/linedet.htm.
%     band=data_toabt.BandGreen;
%     
%     band=ndbi;
 
%% line with a length of three pixels
%     template1=[-1 1 0;
%                -1 1 0;
%                -1 1 0;];
%     template2=[0 1 -1;
%                0 1 -1;
%                0 1 -1;];
%     line_enhanced1 = imfilter(band,template1);
%     line_enhanced2 = imfilter(band,template2);
%     line_enhanced=max(line_enhanced1,line_enhanced2);
%     
%     template1=[-1 -1 -1;
%                 1  1  1;
%                 0  0  0;];
%     template2=[0  0  0;
%                1  1  1;
%               -1 -1 -1;];
%     line_enhanced1 = imfilter(band,template1);
%     line_enhanced2 = imfilter(band,template2);
%     line_enhanced=max(line_enhanced1,line_enhanced);
%     line_enhanced=max(line_enhanced2,line_enhanced);
%     
%     template1=[1 -1 -1;
%                0  1 -1;
%                0  0  1;];
%     template2=[1  0  0;
%               -1  1  0;
%               -1 -1  1;];
%     line_enhanced1 = imfilter(band,template1);
%     line_enhanced2 = imfilter(band,template2);
%     line_enhanced=max(line_enhanced1,line_enhanced);
%     line_enhanced=max(line_enhanced2,line_enhanced);
%     
%     template1=[-1 -1  1;
%                -1  1  0;
%                 1  0  0;];
%     template2=[0  0  1;
%                0  1 -1;
%                1 -1 -1;];
%           
%     line_enhanced1 = imfilter(band,template1);
%     line_enhanced2 = imfilter(band,template2);
%     line_enhanced=max(line_enhanced1,line_enhanced);
%     line_enhanced=max(line_enhanced2,line_enhanced); 
% 	line_enhanced = line_enhanced./3;
    

    template=[-1 2 -1;
              -1 2 -1;
              -1 2 -1;];
	template = template./6;
    line_enhanced = imfilter(band,template);
    
    template=[-1 -1 -1;
               2  2  2;
              -1 -1 -1;];
	template = template./6;
    line_enhanced_new = imfilter(band,template);
    line_enhanced=max(line_enhanced_new,line_enhanced); 
    
    template =[2 -1 -1;
               -1  2 -1;
               -1  -1  2;];
	template = template./6;
    line_enhanced_new = imfilter(band,template);
    line_enhanced=max(line_enhanced_new,line_enhanced); 
 
    template =[-1  -1  2;
               -1  2 -1;
               2 -1 -1;];
	template = template./6;
    line_enhanced_new = imfilter(band,template);
    line_enhanced=max(line_enhanced_new,line_enhanced); 
end

