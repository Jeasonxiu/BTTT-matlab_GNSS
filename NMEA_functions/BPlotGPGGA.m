function [] = BPlotGPGGA(filename)
%
%function [] = BPlotGPGGA(NMEAfile)
%
%   Read the given Logged NMEA file and Plot 'horizonal error', 'Vertical
%   error' from B point
%   
%   input filename : logged NMEA file
%
%   Example : PlotGPGGA('NMEA.txt')
%
%   coded by Joonseong Gim, Jan 7, 2016
%



%% NMEA �α� ���Ͽ��� GPGGA�� gs, la, lo ,h ���� �ϴ� ����
fid=fopen(filename);
%fid=fopen('2016-01-07_113213_B.txt');
check_id = 0; gs_mat = 0; la_mat = 0; lo_mat = 0; h_mat = 0;
while 1
   line = fgetl(fid);
   if ischar(line) == 0, break; end;
   check = char(line(1:6)); % ��ü NMEA �α� ���� �� $GPGGA �� �����Ͽ� ����
    if check == char('$GPGGA')
        gs_temp = str2double(line(7:13));
        gs_mat = [gs_mat; gs_temp]; % ��ü NMEA �α� ������ gs �κ� ����
        la_temp = str2double(line(17:28));
        la_dg = fix(la_temp/100) + (la_temp/100-fix(la_temp/100))*100/60;
        la_mat = [la_mat; la_dg]; % ��ü NMEA �α� ������ la �κ� ����
        lo_temp = str2double(line(31:43));
        lo_dg = fix(lo_temp/100) + (lo_temp/100-fix(lo_temp/100))*100/60;
        lo_mat = [lo_mat; lo_dg]; % ��ü NMEA �α� ������ lo �κ� ����
        H_temp = str2double(line(55:60));
        N_temp = str2double(line(63:67));
        h_temp = H_temp + N_temp;   % Ÿ��ü�� ���
        h_mat = [h_mat; h_temp];    
    end
gs = gs_mat(2:length(gs_mat));
la = la_mat(2:length(la_mat));
lo = lo_mat(2:length(lo_mat));
h = h_mat(2:length(h_mat));
end
gd = [la lo h];

%% gd ����� xyz ��ķ� ��ȯ�ϴ� ����
i = 0;
for i = 1:length(gd)
    xyz(i,:) = gd2xyz(gd(i,:));
    i = i + 1;
end

%% ������ġ ����
real_x = -3041241.741;  % �뼺�������� ���� B ���� x
real_y = 4053944.143;   % �뼺�������� ���� B ���� y
real_z = 3859873.640;   % �뼺�������� ���� B ���� z
real = [real_x real_y real_z];  
real_dg = xyz2gd(real); % �뼺�������� ���� B ���� gd

%% ������ ���Ͽ� x, y, z ���̰� ���
i = 0;
for i = 1:length(xyz)
    xyzs(i,:) = real - xyz(i,:);
    i = i + 1;
end

%% ������ ��
i = 0;
for i = 1:length(xyzs)
    topo(i,:) = xyz2topo(xyzs(i,:), real_dg(1), real_dg(2));
    i = i + 1;
end

%% horizonal error Plot�� max ������ ���� ����
max_north = max(topo(:,1)); % north error�� �ִ밪
max_east = max(topo(:,2));  % east error�� �ִ밪

if max_north >= max_east   % north�� east error ���� ���Ͽ� ū���� ����
    Max = max_north + 1;
    Max = fix(Max); 
else
    Max = max_east + 1;
    Max = fix(Max);
end

%% vertical error Plot�� max ������ ���� ����
max_h = max(topo(:,3))+1; % h error�� �ִ밪
min_h = min(topo(:,3))-1; % h error�� �ּҰ�

%% NMEA �յ� ������ ����
topo_sort = length(topo);
i = 0;
for i = 1:length(topo)  % h < 4 �Ǵ� epoch ����
    if topo(i,3) < 4
        cut = i;
        i = i + 1;
        break;
    end
end
i = 0;
for i = 1:length(topo)-cut  % h < 4 �Ǵ� epoch ���� topo ��� ����
    topo_cut(i,1) = topo(cut+i,1);
    topo_cut(i,2) = topo(cut+i,2);
    topo_cut(i,3) = topo(cut+i,3);
end

%% horizonal error Plot(topo_cut)�� max ������ ���� ����
max_north_cut = max(topo(:,1)); % north error(cut)�� �ִ밪
max_east_cut = max(topo(:,2));  % east error(cut)�� �ִ밪

if max_north_cut >= max_east_cut   % north(cut)�� east error(cut) ���� ���Ͽ� ū���� ����
    Max_cut = max_north_cut + 1;
    Max_cut = fix(Max_cut); 
else
    Max_cut = max_east_cut + 1;
    Max_cut = fix(Max_cut);
end

%% vertical error Plot(topo_cut)�� max ������ ���� ����
max_h = max(topo(:,3))+1; % h error�� �ִ밪
min_h = min(topo(:,3))-1; % h error�� �ּҰ�

% %% NMEA GPGGA ��ü ��� plot 
% figure(101)   
% subplot(2,2,1)  % horizonal error
% plot(topo(:,2),topo(:,1),'b.','markersize',12)
% grid on
% % axis([-Max Max -Max Max])
% axis square
% % set(gca,'XTick',[-Max:1:Max],'YTick',[-Max:1:Max])
% title(['\fontsize{16}horizonal error from B Point '])
% xlabel('\fontsize{12}east(m)')
% ylabel('\fontsize{12}north(m)')    
% 
% subplot(2,2,2)  % vertical error
% % plot(gs(:,1), topo(:,3),'r.:','markersize',10) 
% plot(topo(:,3),'b.:','markersize',10) 
% % axis([0 length(topo) min_h max_h])
% % axis square
% % set(gca,'XTick',[0:200:length(topo)])
% title({'\fontsize{16}vertical error'})
% xlabel(['\fontsize{12}epoch', ' (',num2str(length(topo)),')'])
% ylabel('\fontsize{12}Up(m)')  
% grid on
% 
% % figure(102) % google map plot
% % hold on
% % plot(gd(:,2), gd(:,1),'r.','markersize',10)
% % plot(real_dg(2),real_dg(1),'b.','markersize',13)
% % axis([126.876 126.878 37.479 37.480])
% % axis equal
% % plot_google_map;
% 
% %% NMEA GPGGA �յ� cut Plot
% subplot(2,2,3)  % horizonal error
% plot(topo_cut(:,2),topo_cut(:,1),'b.','markersize',12)
% grid on
% % axis([-Max_cut Max_cut -Max_cut Max_cut])
% axis square
% % set(gca,'XTick',[-Max_cut:1:Max_cut],'YTick',[-Max_cut:1:Max_cut])
% title(['\fontsize{16}horizonal error(cut) from B Point '])
% xlabel('\fontsize{12}east(m)')
% ylabel('\fontsize{12}north(m)')    
% 
% subplot(2,2,4)  % vertical error
% % plot(gs(:,1), topo(:,3),'r.:','markersize',10) 
% plot(topo_cut(:,3),'b.:','markersize',10)
% % axis([0 length(topo_cut) min_h max_h])
% % axis square
% % set(gca,'XTick',[0:200:length(topo_cut)])
% title({'\fontsize{16}vertical error(cut)'})
% xlabel(['\fontsize{12}epoch', ' (',num2str(length(topo_cut)),')'])
% ylabel('\fontsize{12}Up(m)')  
% grid on

% figure(202) % google map plot
% hold on
% plot(gd(:,2), gd(:,1),'r.','markersize',10)
% plot(real_dg(2),real_dg(1),'b.','markersize',13)
% axis([126.876 126.878 37.479 37.480])
% axis equal
% plot_google_map;
%% NMEA GPGGA ��ü ��� plot 
%% NMEA GPGGA ��ü ��� plot 
figure(101)   
subplot(1,2,1)  % horizonal error
hold on
plot(topo(:,2),topo(:,1),'b.','markersize',12)
plot(topo(1,2),topo(1,1),'ro','markersize',7)
plot(topo(length(topo),2),topo(length(topo),1),'go','markersize',7)
legend('S7', 'first', 'last')
grid on
% axis([-Max Max -Max Max])
axis([-10 10 -10 10])
axis square
% set(gca,'XTick',[-Max:1:Max],'YTick',[-Max:1:Max])
title(['\fontsize{16}horizonal error from A Point '])
xlabel('\fontsize{12}east(m)')
ylabel('\fontsize{12}north(m)')    

subplot(1,2,2)  % vertical error
% plot(gs(:,1), topo(:,3),'r.:','markersize',10) 
plot(topo(:,3),'b.:','markersize',10) 
% axis([0 length(topo) min_h max_h])
axis square
% set(gca,'XTick',[0:200:length(topo)])
title({'\fontsize{16}vertical error'})
xlabel(['\fontsize{12}epoch', ' (',num2str(length(topo)),')'])
ylabel('\fontsize{12}Up(m)')  
grid on

% figure(102) % google map plot
% hold on
% plot(gd(:,2), gd(:,1),'r.','markersize',10)
% plot(real_dg(2),real_dg(1),'b.','markersize',13)
% axis([126.876 126.878 37.479 37.480])
% axis equal
% plot_google_map;

