function DATA_STRUCT = rtcm_mode_jprt(fid, f41)
%
%function [] = rtcm_mode_jprt()
%
%   jprt format : { binary data + $GNRMC }<CR>
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
global next_line;   % 다은 데이터 문자열
global now_word;    % 현재 데이터 워드
global sync;        % sync가 맞다면 한 프레임을 빠른속도로 읽음
global start;       % 해당 함수의 첫 호출을 구분하기위함
%% 초기설정 ( SET now_line )
if ~start
    now_line=strtok(fgets(fid), 13);        % <CR> : 13
    start=1;
end
%% JPRT 데이터에서 데이터/시간 분리
next_line = strtok(fgets(fid), 13);         % <CR> : 13
cut = strfind(now_line, '$G');              % 데이터 / 시간 분리
[year month day hour minute second] = ReadRMC(now_line(cut(end):end));
now_line = now_line(1:cut(end)-1);
DATA_STRUCT.time = [year month day hour minute second];
%% 데이터 추출
while 1                         % 한 line을 전부 읽음
    %--- 메세지 첫비트의 위치를 찾음 ----------------------------------------
    while ~sync
        rtcm_find_sync()
    end
    %--- 각 메세지의 해더를 읽어옴 -----------------------------------------
    if length(now_line) > 2
        if rtcm_preamble(now_word)
            %--- 각 메세지의 해더를 읽어옴 ---------------------------------
            header = rtcm_get_header();
            
            %--- type 분리 ------------------------------------------------
            switch header(2)    % message type
                case 1
                    DATA_STRUCT.type1.header = header;
                    DATA_STRUCT.type1.data = rtcm_get_type1(header(6));
%                     bin = rtcm_get_bin(header(6));
%                     val = bin(:, 3:26);
                case 31
                    DATA_STRUCT.type31.header = header;
                    DATA_STRUCT.type31.data = rtcm_get_type31(header(6));
%                     bin = rtcm_get_bin(header(6));
%                     val = bin(:, 3:26);
                case 41
                    rtcm_get_type41(header(6), header(4), DATA_STRUCT.time, f41);
%                     bin = rtcm_get_bin(header(6));
%                     val = bin(:, 3:26);
                otherwise
                    for i=1:header(6) % no. of data wds
                        rtcm_next_word()
                    end
            end
        else                    % sync 잃음 ㅠㅠ
            sync = 0;
        end
    else
        now_line = next_line;
        return
    end
end
    