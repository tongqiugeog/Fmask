function [x_new,y_new]=getRealCloudPositionS2(x,y,h,VZAxy,VAAxy,resolu)
    % imput "x",j col
    % imput "y",i row
    % imput cloud height "h"
    % H=786000; % average S2 height (m)

%     dist_move = h * tan(VZAxy) ./ resolu ; % but the height of cloud's surface should be excluded.
% 
%     % dist=(A*x+B*y+C)/((A^2+B^2)^(0.5));% from the cetral perpendicular (unit: pixel)
%     % dist_par=dist/cos(omiga_per-omiga_par);
%     % dist_move=dist_par.*h/H; % cloud move distance (m);
%     delt_x=dist_move(2)*cos(pi/2-VAAxy);
%     delt_y=dist_move(1)*-sin(pi/2-VAAxy);

    dist_move = h .* tan(VZAxy) ./ resolu ; % but the height of cloud's surface should be excluded.

    delt_x=dist_move(2).*cos(pi/2-VAAxy);
    delt_y=dist_move(1).*-sin(pi/2-VAAxy);

    x_new=x+delt_x; % new x, j
    y_new=y+delt_y; % new y, i
end