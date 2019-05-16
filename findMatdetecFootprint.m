function Matdetec2 = findMatdetecFootprint(DetectFootPrint,XPsizeOut,YPsizeOut)
clear Matdetec k
for i = 1:length(DetectFootPrint.Ndect)
  
  IDX = knnsearch(XPsizeOut',DetectFootPrint.Data{i,1}(:,1));
  IDY = knnsearch(YPsizeOut,DetectFootPrint.Data{i,1}(:,2));
  
  dum2 = single(poly2mask(double(IDX), double(IDY),length(XPsizeOut),length(XPsizeOut))) ;
  clear IDX IDY;
  dum2 = conv2(dum2,ones(3),'same')>0; % to fill boundary
  Matdetec(:,:,i) = dum2;
  clear dum*
  
  %   find orientation of detect + slope for computing perpendicular kernel
  I=nan(size(Matdetec,2),1);
  for ii=1:size(Matdetec,1)
    dum=find(Matdetec(ii,:,i)==1,1,'first');
    if ~isempty(dum)
      I(ii,1)=dum;
    end
  end
  clear dum;
  J = [1:size(Matdetec,1)]' ;
  test = ~isnan(I) & I > 1 & I < size(Matdetec,2);
  warning off all % if warning => not enough point => slope=0 => good because tile boundary
  k{i,1} = polyfit(J(test),I(test),1);
  clear test;
    
  I=nan(size(Matdetec,2),1);
  for ii=1:size(Matdetec,1)
    dum=find(Matdetec(ii,:,i)==1,1,'last');
    if ~isempty(dum)
      I(ii,1)=dum;
    end
  end
  J = [1:size(Matdetec,1)]' ;
  test = ~isnan(I) & I > 1 & I < size(Matdetec,2);
  k{i,2} = polyfit(J(test),I(test),1);
  clear test;
  warning on all

end

% mediane
for i = 1:length(DetectFootPrint.Ndect)-1
  mediane = mean( [k{i,2} ; k{i+1,1}] ) ;
  k{i,2} = mediane ;
  k{i+1,1} = mediane ;
  clear mediane;
end
J = [1:size(Matdetec,1)]' ;
I = [1:size(Matdetec,2)] ;

[Jmat Imat] = meshgrid(I,J);
clear I J;

Matdetec2 = nan(size(Matdetec,1),size(Matdetec,2));
clear Matdetec;
for i = 1:length(DetectFootPrint.Ndect)
  
  liminf = polyval(k{i,1},Jmat);
  limsup = polyval(k{i,2},Jmat);
  
  Matdetec2(Imat>=liminf & Imat<=limsup) = i ;
  clear liminf limsup;
end
clear Imat ImatJ k;

Matdetec2 = Matdetec2';

return

%%
%% old codes
%%

 clear Matdetec k
  for i = 1:length(DetectFootPrint.Ndect)
    
    IDX = knnsearch(XPsizeOut',DetectFootPrint.Data{i,1}(:,1));
    IDY = knnsearch(YPsizeOut,DetectFootPrint.Data{i,1}(:,2));
    
    dum2 = single(poly2mask(double(IDX), double(IDY),length(XPsizeOut),length(XPsizeOut))) ;
    dum2 = conv2(dum2,ones(3),'same')>0; % to fill boundary
    Matdetec(:,:,i) = dum2;
    clear dum*
    
    %   find orientation of detect + slope for computing perpendicular kernel
    I=nan(size(Matdetec,2),1);
    for ii=1:size(Matdetec,1)
      dum=find(Matdetec(ii,:,i)==1,1,'first');
      if ~isempty(dum)
        I(ii,1)=dum;
      end
    end
    J = [1:size(Matdetec,1)]' ;
    test = ~isnan(I) & I > 1 & I < size(Matdetec,2);
    k{i,1} = polyfit(J(test),I(test),1);
    
    I=nan(size(Matdetec,2),1);
    for ii=1:size(Matdetec,1)
      dum=find(Matdetec(ii,:,i)==1,1,'last');
      if ~isempty(dum)
        I(ii,1)=dum;
      end
    end
    J = [1:size(Matdetec,1)]' ;
    test = ~isnan(I) & I > 1 & I < size(Matdetec,2);
    k{i,2} = polyfit(J(test),I(test),1);
    
    
    
    %   n = 100 ;
    %   x =1:n ;
    %   y = round(x*k(1)) ;
    %   y = y-min(y)+1 ;
    %   kern = zeros([max(y) n]) ;
    %   linearInd = sub2ind(size(kern), y, x);
    %   kern(linearInd)=1;
    %   kern = flipud(kern);
    
  end
  
  for i = 1:length(DetectFootPrint.Ndect)-1
    mediane = mean( [k{i,2} ; k{i+1,1}] ) ;
    k{i,2} = mediane ;
    k{i+1,1} = mediane ;
  end
  
  figure
  imagesc(Matdetec(:,:,1))
  hold on
  for i = 1:length(DetectFootPrint.Ndect)
    J = [1:size(Matdetec,1)]' ;
plot(polyval(k{i,2},J),J,'.g')
 plot(polyval(k{i,1},J),J,'r')
  end
  %   find orientation of detect + slope for computing perpendicular kernel
  %   [dum maxdetec]=max(sum(Myreshape(Matdetec,-1)));
  %   I=nan(size(Matdetec,2),1);
  %   for i=1:size(Matdetec,1)
  %     dum=find(Matdetec(i,:,maxdetec)==1,1,'first');
  %     if ~isempty(dum)
  %       I(i,1)=dum;
  %     end
  %   end
  %   J = [1:size(Matdetec,1)]' ;
  %   k = polyfit(J(~isnan(I)),I(~isnan(I)),1);
  %
  %   n = 100 ;
  %   x =1:n ;
  %   y = round(x*k(1)) ;
  %   y = y-min(y)+1 ;
  %   kern = zeros([max(y) n]) ;
  %   linearInd = sub2ind(size(kern), y, x);
  %   kern(linearInd)=1;
  %   kern = flipud(kern);
  
  
  Matdetec = cat(1,Matdetec(1:100,:,:),Matdetec,Matdetec(end-99:end,:,:));
  Matdetec = cat(2,Matdetec(:,1:100,:),Matdetec,Matdetec(:,end-99:end,:));
  
  %   Matdetec=imresize(Matdetec,100/Psize,'bilinear'); %100m => 10m
  %   for i = 1:length(DetectFootPrint.Ndect)
  %   Matdetec2(:,:,i)=conv2(Matdetec(:,:,i),ones(100),'same'); %100m => 10m
  %   end
%   for i = 1:length(DetectFootPrint.Ndect)
%     Matdetec(:,:,i)=conv2(Matdetec(:,:,i),kern,'same'); %100m => 10m
%   end
%   Matdetec = Matdetec(101:end-100,101:end-100,:);
  %   [dum Matdetec2]=max(Matdetec2,[],3);
  %   Matdetec2(dum==0)=NaN;
  [dum Matdetec]=max(Matdetec,[],3);
