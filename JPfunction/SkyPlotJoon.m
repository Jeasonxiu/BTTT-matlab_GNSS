% function SkyPlotJoon(data)
%
%
% input : data(n X 4)
% data(:,1) = prn, data(:,2) = el, data(:,3) = az, data(:,4) = SNR

    clear all; close all;

    load('rover_B_1.mat');

    data = [NMEAQM(:,5:8)];
data = data(data(:,1) > 0, :);
Az = data(:,3); El = data(:,2);
yy = (El-90).* -(cos(Az*pi/180));
xx = (El-90).* -(sin(Az*pi/180));
data(:,2) = xx; data(:,3) = yy;
%     data(:,4) = fix(data(:,4) - 15);
cutoff = 15;

figure(60)
hold on; grid on;
Skymap(cutoff);

if length(data(1,:)) == 4
map = jet(60);
colormap(map)
end

prn = unique(data(:,1))';

for i=1:length(prn)
    prn_ = prn(i);
    
    data_ = data(data(:,1) == prn_, :);
    len = length(data_(:,1));
    if min(prn) < 100
        if prn_ < 10
            PRN = strcat('G0',num2str(prn_));
        elseif prn_ >70
            PRN = strcat('GL',num2str(prn_));
        else
            PRN = strcat('G',num2str(prn_));
        end
    else
        if prn_ < 110
            PRN = strcat('G0',num2str(prn_-100));
        elseif prn_ < 200
            PRN = strcat('G',num2str(prn_-100));
        elseif prn_ < 210
            PRN = strcat('C0',num2str(prn_-200));
        elseif prn_ > 210
            PRN = strcat('C',num2str(prn_-200));
        end
    end
    
    xlast = data_(len,2);
    ylast = data_(len,3);
    if length(data(1,:)) == 4
        for m = 1:len
            color = map(data_(m,4),:);
            x_ = data_(m,2); y_ = data_(m,3);
            plot(x_, y_, 'o','color', color);
            plot(x_, y_, '*','color', color);
            text(xlast, ylast, PRN,'Fontsize',15)
        end
    elseif length(data(1,:)) == 3
        x_ = data_(:,2); y_ = data_(:,3);
        if prn_ < 206
            plot(x_, y_, 'r.','Markersize',20);
            text(xlast, ylast, PRN,'Fontsize',15)
        elseif prn_ == 206 | prn_ == 207 | prn_ == 208 | prn_ == 209 | prn_ == 210 | prn_ == 213
            plot(x_, y_, 'g.','Markersize',20);
            text(xlast, ylast, PRN,'Fontsize',15)
        else
            plot(x_, y_, 'b.','Markersize',20);
            text(xlast, ylast, PRN,'Fontsize',15)
        end
    end
end


if length(data(1,:)) == 4
colorbar('Ticks',[0.25, 0.5, 0.75],...
    'Ticklabels',{'30', '40', '50'});
end




