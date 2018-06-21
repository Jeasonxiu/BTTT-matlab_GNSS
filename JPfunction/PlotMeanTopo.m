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

%% user_mean ����
x_mean = mean(user_xyz(:,1)); y_mean = mean(user_xyz(:,2)); z_mean = mean(user_xyz(:,3));   % user�� x,y,z �� ������ ���
user_mean = [x_mean y_mean z_mean];
%% ��հ��� ���Ͽ� x, y, z ���̰� ���
for j = 1:length(user_xyz)
    xyzs(j,:) = user_mean - user_xyz(j,:);
end
user_mean_dg = xyz2gd(user_mean);
%% ���̰��� �̿��Ͽ� delta_n, delta_e, h���
for k = 1:length(xyzs)
    topo(k,:) = xyz2topo(xyzs(k,:), user_mean_dg(1), user_mean_dg(2));
end

%% horizonal error Plot�� max ������ ���� ����
max_north = abs(max(topo(:,1))); % north error�� �ִ밪
max_east = abs(max(topo(:,2)));  % east error�� �ִ밪

if max_north >= max_east   % north�� east error ���� ���Ͽ� ū���� ����
    Max = max_north + 1;
    Max = fix(Max); 
else
    Max = max_east + 1;
    Max = fix(Max);
end

if Max > 10
    Max = round(Max/10)*10;
end
%% xtick ���� ����
if Max <= 10
    tick = 1;
elseif Max <= 20
    tick = 2;
elseif Max <= 30
    tick = 3;
else
    tick = 5;
end


%% vertical error Plot�� max ������ ���� ����
max_h = max(topo(:,3))+1; % h error�� �ִ밪
min_h = min(topo(:,3))-1; % h error�� �ּҰ�

%% user_mean �� ���� Topology Plot
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