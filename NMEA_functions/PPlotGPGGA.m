function [] = PPlotGPGGA(filename)
%
%function [] = PPlotGPGGA(NMEAfile)
%
%   Read the given Logged NMEA file and Plot 'horizonal error', 'Vertical
%   error', from selecitd point
%   
%   input filename : logged NMEA file
%
%   Example : PPlotGPGGA('NMEA.txt','a')
%
%   coded by Joonseong Gim, Jan 7, 2016
%



%% NMEA 로그 파일에서 GPGGA의 gs, la, lo ,h 추출 하는 과정
% point = 'b';
% filename = 'DBUU049i_16.ubx';
point = filename(4);
% point = 'B';
fid=fopen(filename,'r');
if fid == -1
    disp('Cannot locate the input file!')
else
    GPGGA = getGPGGA(filename);
    GPRMC = getGPRMC(filename);
    name = filename(length(filename)-5:length(filename)-4);
    [gd, gs, utc, la, lo, h, xyz] = GGA2gd(filename);
    if ~isempty(GPRMC)
        for ii = 1: length(GPRMC)
            line = cell2mat(GPRMC(ii,1));
            if length(line) >= 70
                break
            end
        end
        if length(line) >= 60
            index = findstr(line,',');
            yymmdd = line(index(9)+1:index(10)-1);
            if str2num(yymmdd(5:6)) >= 80
                yyyy = 1900 + str2num(yymmdd(5:6));
            else
                yyyy = 2000 + str2num(yymmdd(5:6));
            end
            mm = str2num(yymmdd(3:4));
            dd = str2num(yymmdd(5:6));
        else
            yyyy = 0;
            mm = 0;
            dd = 0;
        end
    else
        yyyy = 0;
        mm = 0;
        dd = 0;
    end
    
    %% 선택 기준위치 참값
    [real, real_dg] = choiceP(point);
%     cut = round(length(xyz)/10);
    if gs(1,1) == 0
        [dXYZ, dNEV] = PosTErrorsGPGGA(utc(:,1), real, xyz,name,yyyy,mm,dd,point);
%         [dXYZ, dNEV] = PosTErrorsGPGGA(utc(cut:length(xyz)-cut,1), real, xyz(cut:length(xyz)-cut,1:3),name,yyyy,mm,dd,point);
    else
        [dXYZ, dNEV] = PosTErrorsGPGGA(gs(:,1), real, xyz,name,yyyy,mm,dd,point);
%         [dXYZ, dNEV] = PosTErrorsGPGGA(gs(cut:length(xyz)-cut,1), real, xyz(cut:length(xyz)-cut,1:3),name,yyyy,mm,dd,point);
    end
    
end




% %% 참값과 비교하여 x, y, z 차이값 계산
% i = 0;
% for i = 1:length(xyz)
%     xyzs(i,:) = real - xyz(i,:);
%     i = i + 1;
% end
% 
% %% 참값과 비교하여 delta_n, delta_e, h계산
% i = 0;
% for i = 1:length(xyzs)
%     topo(i,:) = xyz2topo(xyzs(i,:), real_dg(1), real_dg(2));
%     i = i + 1;
% end
% 
% %% horizonal error Plot의 max 선정을 위한 과정
% [pick,  value, direction] = picktopo(topo)
% 
% %% NMEA 앞뒤 데이터 삭제
% cut = fix(length(topo)/10);
% i = 0;
% for i = 1:length(topo)-2*cut  % h < 4 되는 epoch 부터 topo 행렬 생성
%     topo_cut(i,1) = topo(cut+i,1);
%     topo_cut(i,2) = topo(cut+i,2);
%     topo_cut(i,3) = topo(cut+i,3);
%     gd_cut(i,:) = gd(cut+i,:);
% end
% 
% %% NMEA GPGGA 전체 결과 plot 
% figure(101)   
% subplot(1,2,1)  % horizonal error
% plot(topo(:,2),topo(:,1),'ro','MarkerEdgeColor','b','MarkerFaceColor','r','markersize',4)
% grid on
% axis square
% axis ([-pick pick -pick pick])
% title(['\fontsize{16}horizonal error from Point ',point])
% xlabel('\fontsize{12}east(m)')
% ylabel('\fontsize{12}north(m)')    
% 
% subplot(1,2,2)  % vertical error
% plot(topo(:,3),'r.:','markersize',5) 
% axis square
% xlim([0 length(topo)])
% set(gca,'XTick',[0:200:length(topo)])
% title({'\fontsize{16}vertical error'})
% xlabel(['\fontsize{12}epoch', ' (',num2str(length(topo)),')',' 1Hz'])
% ylabel('\fontsize{12}Up(m)')  
% grid on
% 
% % figure(201) % google map plot
% % hold on
% % plot(gd(:,2), gd(:,1),'r.','markersize',10)
% % plot(real_dg(2),real_dg(1),'b.','markersize',13)
% % axis([126.876 126.878 37.479 37.480])
% % axis equal
% % plot_google_map;
% 
% %% NMEA GPGGA 앞뒤 cut Plot
% figure(102)
% subplot(1,2,1)  % horizonal error
% plot(topo_cut(:,2),topo_cut(:,1),'ro','MarkerEdgeColor','b','MarkerFaceColor','r','markersize',4)
% grid on
% axis square
% title(['\fontsize{16}horizonal error(cut) from Point ',point])
% xlabel('\fontsize{12}east(m)')
% ylabel('\fontsize{12}north(m)')    
% 
% subplot(1,2,2)  % vertical error
% plot(topo_cut(:,3),'r.:','markersize',5)
% xlim([0 length(topo_cut)])
% axis square
% set(gca,'XTick',[0:200:length(topo_cut)])
% title({'\fontsize{16}vertical error(cut)'})
% xlabel(['\fontsize{12}epoch', ' (',num2str(length(topo_cut)),')',' 1Hz'])
% ylabel('\fontsize{12}Up(m)')  
% grid on

% figure(202) % google map plot
% hold on
% plot(gd(:,2), gd(:,1),'r.','markersize',10)
% plot(real_dg(2),real_dg(1),'b.','markersize',13)
% axis([126.876 126.878 37.479 37.480])
% axis equal
% plot_google_map;

%% gd의 latitude, longitude plot
% figure(103)
% subplot(4,1,1)  % latitude plot whole data
% plot(gd(:,1),'r.:','markersize',8) 
% xlim([0 length(gd)])
% set(gca,'XTick',[0:200:length(gd)])
% title({'\fontsize{16}Latitude'})
% xlabel(['\fontsize{12}epoch', ' (',num2str(length(gd)),')',' 1Hz'])
% ylabel('\fontsize{12}Degree(°)')  
% grid on 
% subplot(4,1,2)  % latitude plot whole data
% plot(gd(:,2),'r.:','markersize',8) 
% xlim([0 length(gd)])
% set(gca,'XTick',[0:200:length(gd)])
% title({'\fontsize{16}Longitude'})
% xlabel(['\fontsize{12}epoch', ' (',num2str(length(gd)),')',' 1Hz'])
% ylabel('\fontsize{12}Degree(°)')  
% grid on 
% subplot(4,1,3)  % latitude plot whole data
% plot(gd_cut(:,1),'b.:','markersize',8) 
% xlim([0 length(gd_cut)])
% set(gca,'XTick',[0:200:length(gd_cut)])
% title({'\fontsize{16}Latitude'})
% xlabel(['\fontsize{12}epoch', ' (',num2str(length(gd_cut)),')',' 1Hz'])
% ylabel('\fontsize{12}Degree(°)')  
% grid on 
% subplot(4,1,4)  % latitude plot whole data
% plot(gd_cut(:,2),'b.:','markersize',8) 
% xlim([0 length(gd_cut)])
% set(gca,'XTick',[0:200:length(gd_cut)])
% title({'\fontsize{16}Longitude'})
% xlabel(['\fontsize{12}epoch', ' (',num2str(length(gd_cut)),')',' 1Hz'])
% ylabel('\fontsize{12}Degree(°)')  
% grid on 


