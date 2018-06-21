function DATA_STRUCT = rtcm_mode_ntrip(fid)
%
%function [] = rtcm_mode_ntrip()
%
%   ntrip format : { time[yyyy-mm-ddThh:mm:ss]<CR> + binary data }<CR>
%
%   출력 데이터 구조는 아래와같고, 각 필드의 데이터 순서는 RTCM 문서를 따름
%
%   <output>
%   DATA_STRUCT -time
%               -type## -header
%                       -data
%               ...
%
%   Copyright: taeil Kim, February 20, 2015@INHA University

%% Global
global now_line;    % 현재 데이터 문자열
global next_line;   % 다음 데이터 문자열
global now_tline;   % 현재 시간 문자열
global next_tline;  % 다음 시간 문자열
global now_word;    % 현재 데이터 워드
global sync;        % sync가 맞다면 한 프레임을 빠른속도로 읽음
global start;       % 해당 함수의 첫 호출을 구분하기위함
global glo_hour; 
global glo_min; 
global glo_sec; 
%% 초기설정 ( SET now_line )
if ~start
    now_tline=fgets(fid);
    while length(now_tline) < 18    % 첫번째 줄이 공백일때의 예외처리
        now_tline = fgets(fid);
    end
    now_line =strtok(fgets(fid), 13);       % <CR> : 13
    start=1;
end
%% NTRIP 데이터에서 데이터/시간 분리
next_tline= fgets(fid);
next_line = strtok(fgets(fid), 13);         % <CR> : 13
year  = str2num( now_tline( 1: 4) );
month = str2num( now_tline( 6: 7) );
day   = str2num( now_tline( 9:10) );
hour  = str2num( now_tline(12:13) );
minute= str2num( now_tline(15:16) );
second= str2num( now_tline(18:19) );
DATA_STRUCT.time = [year month day hour minute second];
glo_hour = hour;
glo_min = minute;
glo_sec = second;
%% 데이터 추출
while 1                         % 한 line을 전부 읽음
    %--- 메세지 첫비트의 위치를 찾음 ----------------------------------------
    while ~sync
        rtcm_find_sync()
    end
    %--- 프레임의 첫 워드가 맞는지 검사 -------------------------------------
    if length(now_line) > 2
        if rtcm_preamble(now_word)
            %--- 각 메세지의 해더를 읽어옴 ---------------------------------
            header = rtcm_get_header();
            
            %--- type 분리 ------------------------------------------------
            switch header(2)    % message type
                case 1
                    DATA_STRUCT.type1.header = header;
                    DATA_STRUCT.type1.data = rtcm_get_type1(header(6));
                case 31
                    DATA_STRUCT.type31.header = header;
                    DATA_STRUCT.type31.data = rtcm_get_type31(header(6));
                otherwise
                    for i=1:header(6) % no. of data wds
                        rtcm_next_word()
                    end
            end
        else                    % sync 잃음 ㅠㅠ
            sync = 0;
        end
    else
        now_tline= next_tline;
        now_line = next_line;
        return
    end
end
    