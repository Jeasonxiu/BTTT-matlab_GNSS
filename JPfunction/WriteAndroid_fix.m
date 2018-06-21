% function [estm_fix] = WriteAndroid_fix(filename)
%
%   function [LLH] = WriteAndroid_fix(filename)
%   Read the given logged Android file and make a 'LLH'
%
%   ex) [LLH] = WriteAndroid_fix'pseudoranges_log_2017_02_14_18_58_59.txt');
%
%   Coded by Joonseong Gim, Sep 18, 2017
%

clear all
close all
% filename = 'pseudoranges_log_2017_02_14_18_58_59.txt';
% filename = 'pseudoranges_log_2017_03_02_15_33_59.txt';      % �������� �޾��ֽ� �ڷ�
% filename = 'pseudoranges_log_2017_08_01_15_39_28.txt';      % �뼺 B ���� ����
% filename = 'pseudoranges_log_2017_08_01_15_50_15.txt';      % �뼺 B ���� �����
% filename = 'pseudoranges_log_2017_08_01_16_02_35.txt';      % �뼺 B ���� ���� ����
% filename = 'pseudoranges_log_2017_08_01_16_14_00.txt';      % �뼺 B ���� ���� ����
% filename = 'pseudoranges_log_2017_09_14_16_31_57.txt';      % �뼺 B ���� 1
% filename = 'pseudoranges_log_2017_09_14_16_42_37.txt';      % �뼺 B ���� 2
% filename = 'pseudoranges_log_2017_09_14_16_53_37.txt';      % �뼺 B ���� 3
% filename = 'pseudoranges_log_2017_09_14_17_04_08.txt';      % �뼺 B ���� 4
% filename = 'pseudoranges_log_2017_09_14_17_14_31.txt';      % �뼺 B ���� 5
% filename = 'pseudoranges_log_2017_09_14_17_25_10.txt';      % �뼺 B ���� 6
% filename = 'pseudoranges_log_2017_09_14_17_36_04.txt';      % �뼺 B ���� 7
% filename = 'pseudoranges_log_2017_09_14_17_46_35.txt';      % �뼺 B ���� 8
% filename = 'pseudoranges_log_2017_09_14_17_57_01.txt';      % �뼺 B ���� 9
% filename = 'pseudoranges_log_2017_09_14_18_08_57.txt';      % �뼺 B ���� 10
% filename = 'pseudoranges_log_2017_09_19_09_39_30.txt';      % �뼺 Round
% filename = 'pseudoranges_log_2017_09_19_09_41_27.txt';      % �뼺 Round
filename = 'pseudoranges_log_2017_10_13_14_46_27.txt';

%% ���� ���
CCC = 299792458;                                    % Speed of Light
WeekSecond = 604800;                                % Week Second
%% Type ���� ���
C1 = 20;                                                                 % Code
C1PrSigmaM = 21;                                                         % PrSigmaM
C1prrSigmaMps = 22;                                                      % PseudorangeRateUncertaintyMetersPerSecond
L1 = 11;                                                                 % reserved
S1 = 41;                                                                 % SNR
U1 = 51;                                                                 % Uncertainty

%% text ���� ���� �� QMfile ����
fid = fopen(filename,'r');

%% header ���� ����
ready = 0;
while ready == 0
    s = fgets(fid);
    if length(s) > 100
        if ~isempty(strfind(s,'ElapsedRealtimeMillis'))
            ready = 1;
        end
    end
end
%% QMfile ���� ���� ����
cnt = 0; cnt2=1; Epoch = 0; Line = 1;
while 1
    line = fgets(fid);
    if length(line) >= 72
        if line(1:3) == 'Fix'
            cnt = cnt + 1;
            fix = strsplit(line,',');       % GGA cell array
            Fix = textscan(line(1:end),'%s%s%f%f%f%f%f%f','delimiter',',');
            formatOut = 'yyyy-mm-dd/HH:MM:SS.FFF';
            date = datestr((cell2mat(Fix(8))+32400000)/86400/1000 + datenum(1970,1,1),formatOut);
            yyyy = str2num(date(1:4));
            mo = str2num(date(6:7));
            dd = str2num(date(9:10));
            hh = str2num(date(12:13));
            mm = str2num(date(15:16));
            ss = str2num(date(18:23));
            [GW, gs] = date2gwgs(yyyy,mo,dd,hh,mm,ss);
            Lat = cell2mat(Fix(3));
            Lon = cell2mat(Fix(4));
            Height = cell2mat(Fix(5));
            xyz = gd2xyz([Lat,Lon,Height]);
            LLH(cnt,:) = [gs, xyz, Lat, Lon, Height];
        end
    end
    %% While �� ����
    if line == -1, break;  end
end
estm_fix = LLH;

plot(estm_fix(:,6), estm_fix(:,5),'bo')
plot_google_map