function [user_mean, user_mean_dg] = PlotMeanTopo(user_xyz)
%
%function [user_mean] = PlotMeanTopo(estm)
%
%   Load the estimated estm and Plot 'horizonal error', 'Vertical
%   error', from user_xyz position(x,y,z)
%   
%   input filename : user_xyz(x y z)(ecef)
%   
%   Example : PlotMeanTopo(user_xyz)
%
%   coded by Joonseong Gim, Jan 13, 2016
%

%% user_mean 생성
x_mean = mean(user_xyz(:,1)); y_mean = mean(user_xyz(:,2)); z_mean = mean(user_xyz(:,3));   % user의 x,y,z 값 각각의 평균
user_mean = [x_mean y_mean z_mean];
%% 평균값과 비교하여 x, y, z 차이값 계산
for j = 1:length(user_xyz)
    xyzs(j,:) = user_mean - user_xyz(j,:);
end
user_mean_dg = xyz2gd(user_mean);
%% 차이값을 이용하여 delta_n, delta_e, h계산
for k = 1:length(xyzs)
    topo(k,:) = xyz2topo(xyzs(k,:), user_mean_dg(1), user_mean_dg(2));
end

%% horizonal error Plot의 max 선정을 위한 과정
max_north = abs(max(topo(:,1))); % north error의 최대값
max_east = abs(max(topo(:,2)));  % east error의 최대값

if max_north >= max_east   % north와 east error 값을 비교하여 큰값을 선택
    Max = max_north + 1;
    Max = fix(Max); 
else
    Max = max_east + 1;
    Max = fix(Max);
end

if Max > 10
    Max = round(Max/10)*10;
end
%% xtick 선정 과정
if Max <= 10
    tick = 1;
elseif Max <= 20
    tick = 2;
elseif Max <= 30
    tick = 3;
else
    tick = 5;
end


%% vertical error Plot의 max 선정을 위한 과정
max_h = max(topo(:,3))+1; % h error의 최대값
min_h = min(topo(:,3))-1; % h error의 최소값

%% user_mean 값 기준 Topology Plot
figure(300)
subplot(1,2,1)  % horizonal error
plot(topo(:,2),topo(:,1),'ro','MarkerEdgeColor','y','MarkerFaceColor','r','markersize',3)
grid on
axis([-Max Max -Max Max])
axis square
set(gca,'XTick',[-Max:tick:Max],'YTick',[-Max:tick:Max])
title(['\fontsize{16}horizonal error from User mean'])
xlabel('\fontsize{12}east(m)')
ylabel('\fontsize{12}north(m)')    

subplot(1,2,2)  % vertical error
plot(topo(:,3),'r.:','markersize',8) 
axis([0 length(topo) min_h max_h])
axis square
set(gca,'XTick',[0:200:length(topo)])
title({'\fontsize{16}vertical error'})
xlabel(['\fontsize{12}epoch', ' (',num2str(length(topo)),')',' 1Hz'])
ylabel('\fontsize{12}Up(m)')  
grid on