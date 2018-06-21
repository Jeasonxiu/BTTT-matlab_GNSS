% function plot_source2dis(source, distance)

% input
%       source : source [N x 3] matrix
%               N is number of source (x,y,z(ECEF))
%       distance : distance from source coordinate(km)
%
%       example : plot_source2dis(source, 260);

source = [-3003051.386, 4059906.022, 3883100.080;
    -3113998.569, 3920042.804, 3938650.149;
    -3177399.906, 4290816.111, 3477703.760;
    -3290252.388, 4018635.966, 3690049.353]; % °­È­ °í¼º, ¼­±ÍÆ÷ ÁÂÇ¥

distance = 260;
angulardis = distance/6373;
countofsource = length(source(:,1));
cnt=1;



for j = 1:4
    gd1 = xyz2gd(source(j,:));

    for i = 1:360;
        deg = i;
        des_la = asind(sind(gd1(1))*cos(angulardis)+cosd(gd1(1))*sin(angulardis)*cosd(deg));
        des_lo = gd1(2) + atan2d(sind(deg)*sin(angulardis)*cosd(gd1(1)), cos(angulardis)-sind(gd1(1))*sind(des_la));
        des_gd(i,1:2) = [des_la, des_lo];
    end
    figure(1)
    plot(gd1(2),gd1(1),'*r','Markersize',20)
    hold on
    plot(des_gd(:,2), des_gd(:,1),'-','Linewidth',1)
    
end

plot_google_map;