function [] = WritePRC(filename, mode)
%
%function [] = WritePRC(filename, mode)
%
%   RTCM ������ binary �����͸� �а�,
%   �Ʒ� �������� DGPS, DGLONASS �����͸� ���
%
%   GPS Week Second / PRN + type / PRC / RRC ( type : GPS-0, GLO-64 )
%       <���� type(Glonass�� 65���� ����)�� NMEA���� ��µǴ� ������ ��������>
%
%   <input>
%       filename    : ���ϸ�
%       mode        : {jprt} | ntrip
%                       jprt  - �������̽� ���� ������
%                       ntrip - ntrip ������
%
%   Copyright: taeil Kim, February 20, 2015@INHA University

%% �Է� ���� �˻�
error(nargchk(1, 2, nargin))
if nargin < 2, mode = 'jprt'; end
%% ������ ����
fprintf(1,'������ ������...\n');
fid = fopen('PRCfile', 'w');
DATA = getRTCM(filename, mode);
len = length(DATA);
%% File open

%% ������ ���
fprintf(1,'������ �����...\n');
for i=1:len
    time = DATA{i}.time;
    %--- �ð� ���ذ� ��ȯ --------------------------------------------------
    switch mode
        case 'jprt'         % UTC -> GPST
            [year month day hour minute second] = utc2gpst(time);
            [dum time] = date2gwgs(year, month, day, hour, minute, second);
        case 'ntrip'        % KST -> GPST
            [year month day hour minute second] = kst2utc(time);
            [year month day hour minute second] = ...
                utc2gpst([year month day hour minute second]);
            [dum time] = date2gwgs(year, month, day, hour, minute, second);
        otherwise
            break;
    end
    %--- DGPS Corrections -------------------------------------------------
    if sum( strcmp(fieldnames(DATA{i}), 'type1') )
        data = sortrows( DATA{i}.type1.data, 3 );   % prn ������ ����
        d_len = length(data(:,1));
        s_time= rtcm_time_corr(time, DATA{i}.type1.header(4), 'GPS');
        for j=1:d_len
            % GPS Week Second / prn / PRC / RRC
            fprintf(fid,'%8.1f %5d %10.2f %8.3f \n', s_time, data(j,3), ...
                bitcmp2( data(j,4), 16 ) * ( 0.02  + data(j,1)*0.3 ), ...
                bitcmp2( data(j,5), 8 ) * ( 0.002 + data(j,1)*0.03 ));
        end
    end
    %--- DGLONASS Corrections ---------------------------------------------
    if sum( strcmp(fieldnames(DATA{i}), 'type31') )
        data = sortrows( DATA{i}.type31.data, 3 );  % prn ������ ����
        d_len = length(data(:,1));
        s_time= rtcm_time_corr(time, DATA{i}.type31.header(4), 'GLONASS');
        for j=1:d_len
            % GPS Week Second / prn + 64 / PRC / RRC
            fprintf(fid,'%8.1f %5d %10.2f %8.3f \n', s_time, data(j,3)+64, ...
                bitcmp2( data(j,4), 16 ) * ( 0.02  + data(j,1)*0.3 ), ...
                bitcmp2( data(j,5), 8 ) * ( 0.002 + data(j,1)*0.03 ));
        end
    end
end
%% ������
fprintf(1,'�Ϸ�\n');
fclose(fid);