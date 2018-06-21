function [] = PlotGPGGA(filename)
%
%function [] = PlotGPGGA(NMEAfile)
%
%   Read the given Logged NMEA file and Plot 'horizonal error', 'Vertical
%   error'
%   
%   input filename : logged NMEA file
%
%   Example : PlotGPGGA('NMEA.txt')
%
%   coded by Joonseong Gim, Jan 7, 2016
%



%% NMEA �α� ���Ͽ��� GPGGA�� gs, la, lo ,h ���� �ϴ� ����
% filename = 'test.txt';
fid=fopen(filename,'r');
if fid == -1
    disp('Cannot locate the input file!')
else
    GPGGA = getGPGGA(filename);
    GPRMC = getGPRMC(filename);
    [gd, gs, utc, la, lo, h, xyz] = GGA2gd(filename);
end
%% user_position mean�� ���� topology ��� ����
    % [user_mean,user_mean_dg] = PlotMeanTopo(xyz);

    figure(102) % google map plot
    hold on
    plot(gd(:,2), gd(:,1),'ro','MarkerEdgeColor','y','MarkerFaceColor','r','markersize',3)
    % plot(user_mean_dg(2), user_mean_dg(1),'b*','MarkerEdgeColor','g','MarkerFaceColor','b','markersize',5)
    axis([(min(gd(:,2))-0.001) (max(gd(:,2))+0.001) (min(gd(:,1))-0.001) (max(gd(:,1))+0.001)])
    axis equal
    plot_google_map;


