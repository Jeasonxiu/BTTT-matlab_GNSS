function DATA_STRUCT = rtcm_mode_jprt(fid, f41)
%
%function [] = rtcm_mode_jprt()
%
%   jprt format : { binary data + $GNRMC }<CR>
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
global now_word;    % ���� ������ ����
global sync;        % sync�� �´ٸ� �� �������� �����ӵ��� ����
global start;       % �ش� �Լ��� ù ȣ���� �����ϱ�����
%% �ʱ⼳�� ( SET now_line )
if ~start
    now_line=strtok(fgets(fid), 13);        % <CR> : 13
    start=1;
end
%% JPRT �����Ϳ��� ������/�ð� �и�
next_line = strtok(fgets(fid), 13);         % <CR> : 13
cut = strfind(now_line, '$G');              % ������ / �ð� �и�
[year month day hour minute second] = ReadRMC(now_line(cut(end):end));
now_line = now_line(1:cut(end)-1);
DATA_STRUCT.time = [year month day hour minute second];
%% ������ ����
while 1                         % �� line�� ���� ����
    %--- �޼��� ù��Ʈ�� ��ġ�� ã�� ----------------------------------------
    while ~sync
        rtcm_find_sync()
    end
    %--- �� �޼����� �ش��� �о�� -----------------------------------------
    if length(now_line) > 2
        if rtcm_preamble(now_word)
            %--- �� �޼����� �ش��� �о�� ---------------------------------
            header = rtcm_get_header();
            
            %--- type �и� ------------------------------------------------
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
        else                    % sync ���� �Ф�
            sync = 0;
        end
    else
        now_line = next_line;
        return
    end
end
    