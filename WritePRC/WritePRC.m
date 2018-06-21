function [] = WritePRC(filename, mode)
%
%function [] = WritePRC(filename, mode)
%
%   RTCM 포멧의 binary 데이터를 읽고,
%   아래 형식으로 DGPS, DGLONASS 데이터를 기록
%
%   GPS Week Second / PRN + type / PRC / RRC ( type : GPS-0, GLO-64 )
%       <위의 type(Glonass는 65부터 시작)은 NMEA에서 출력되는 형식을 따른것임>
%
%   <input>
%       filename    : 파일명
%       mode        : {jprt} | ntrip
%                       jprt  - 지평스페이스 옥상 데이터
%                       ntrip - ntrip 데이터
%
%   Copyright: taeil Kim, February 20, 2015@INHA University

%% 입력 변수 검사
error(nargchk(1, 2, nargin))
if nargin < 2, mode = 'jprt'; end
%% 데이터 추출
fprintf(1,'데이터 추출중...\n');
fid = fopen('PRCfile', 'w');
DATA = getRTCM(filename, mode);
len = length(DATA);
%% File open

%% 데이터 기록
fprintf(1,'데이터 기록중...\n');
for i=1:len
    time = DATA{i}.time;
    %--- 시간 기준계 변환 --------------------------------------------------
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
        data = sortrows( DATA{i}.type1.data, 3 );   % prn 순서로 정렬
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
        data = sortrows( DATA{i}.type31.data, 3 );  % prn 순서로 정렬
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
%% 마무리
fprintf(1,'완료\n');
fclose(fid);