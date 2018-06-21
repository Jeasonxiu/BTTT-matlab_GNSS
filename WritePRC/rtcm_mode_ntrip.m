function DATA_STRUCT = rtcm_mode_ntrip(fid)
%
%function [] = rtcm_mode_ntrip()
%
%   ntrip format : { time[yyyy-mm-ddThh:mm:ss]<CR> + binary data }<CR>
%
%   ��� ������ ������ �Ʒ��Ͱ���, �� �ʵ��� ������ ������ RTCM ������ ����
%
%   <output>
%   DATA_STRUCT -time
%               -type## -header
%                       -data
%               ...
%
%   Copyright: taeil Kim, February 20, 2015@INHA University

%% Global
global now_line;    % ���� ������ ���ڿ�
global next_line;   % ���� ������ ���ڿ�
global now_tline;   % ���� �ð� ���ڿ�
global next_tline;  % ���� �ð� ���ڿ�
global now_word;    % ���� ������ ����
global sync;        % sync�� �´ٸ� �� �������� �����ӵ��� ����
global start;       % �ش� �Լ��� ù ȣ���� �����ϱ�����
global glo_hour; 
global glo_min; 
global glo_sec; 
%% �ʱ⼳�� ( SET now_line )
if ~start
    now_tline=fgets(fid);
    while length(now_tline) < 18    % ù��° ���� �����϶��� ����ó��
        now_tline = fgets(fid);
    end
    now_line =strtok(fgets(fid), 13);       % <CR> : 13
    start=1;
end
%% NTRIP �����Ϳ��� ������/�ð� �и�
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
%% ������ ����
while 1                         % �� line�� ���� ����
    %--- �޼��� ù��Ʈ�� ��ġ�� ã�� ----------------------------------------
    while ~sync
        rtcm_find_sync()
    end
    %--- �������� ù ���尡 �´��� �˻� -------------------------------------
    if length(now_line) > 2
        if rtcm_preamble(now_word)
            %--- �� �޼����� �ش��� �о�� ---------------------------------
            header = rtcm_get_header();
            
            %--- type �и� ------------------------------------------------
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
        else                    % sync ���� �Ф�
            sync = 0;
        end
    else
        now_tline= next_tline;
        now_line = next_line;
        return
    end
end
    