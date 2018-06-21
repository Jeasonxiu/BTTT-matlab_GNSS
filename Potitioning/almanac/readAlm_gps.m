function [alm] = readAlm_gps(filename)
%
% function [alm] = readAlm_gps(filename)
%
% Read the GPS Almanac text file, make a alm matrix
%
% input : GPS almanac text
%
% output : GPS almanac Matrix
%
% Example : [alm] = readAlm_gps('current.txt');
%
% coded by Joonseong Gim, DEC 28, 2016

% clear all; close all;

%% ���� ��¥�� �������� navcen���� �˸��� ������ �ٿ�޾� yyyymmdd_almanac.txt ���·� ����
current_time = clock;
if current_time(2) > 9 
    if current_time(3) > 9
        almanacfilename = strcat(num2str(current_time(1)),num2str(current_time(2)),num2str(current_time(3)),'_almanac.txt');
    else
        almanacfilename = strcat(num2str(current_time(1)),num2str(current_time(2)),'0',num2str(current_time(3))','_almanac.txt');
    end
elseif current_time(3) > 9
    almanacfilename = strcat(num2str(current_time(1)),'0',num2str(current_time(2)),num2str(current_time(3)),'_almanac.txt');
else
    almanacfilename = strcat(num2str(current_time(1)),'0',num2str(current_time(2)),'0',num2str(current_time(3)),'_almanac.txt');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filename = '2016_875.txt' ;       % �Լ� �׽�Ʈ�� ����
% filename = almanacfilename;

%% �˸��� ���� ���� ���� �Ǵ�
fid = fopen(filename,'r');
if fid == -1
    disp('No GPS almanac file!')
    urlwrite('http://www.navcen.uscg.gov/?pageName=currentAlmanac&format=yuma-txt',almanacfilename);
    fid = fopen(filename,'r');
    alm = [];
else
    ready = 0; i=1;j=1;
    while ready == 0
        line=fgetl(fid);
        if line == -1
            ready=1;
        end
        if length(line) > 1 & line(1:3) ~= '***'
            [para, value] = strread(line(1,:),'%s%s','delimiter',':');
            alm(j,i) = str2num(cell2mat(value));
            i=i+1;
            if i == 14
                j=j+1;
                i=1;
            end
            
        end
    end
end
alm(:,13) = alm(:, 13) + 1024;
%% �˸��� �����͸� �̿��� ���� ��ǥ ����ϴ� �Լ�
% [SatPos] = GetSatPos_almanac(alm, 5, 100000);


