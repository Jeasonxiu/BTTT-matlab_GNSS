function [] = rtcm_next_word()
%
%function [] = rtcm_next_word()
%
%   now_word�� ���� ������ �������� �ٲ�
%
%   Copyright: taeil Kim, February 20, 2015@INHA University

%% Global
global now_word;
global now_line;
global next_line;
%%
%--- �� ������ �����Ͱ� �����ٿ� ������ -------------------------------------
len = 5 - length(now_line);

if len > 0
    now_line = [now_line next_line(1:len)];
    next_line(1:len) = [];
end

%--- 5���� ����(data�� 30bit) shift ----------------------------------------
for i=1:5
    u_char = uint32(now_line(i));
    for j=1:6
        now_word = bitshift(now_word, 1) + bitget(u_char, j);
    end
end
%--- ����� ���� ���� --------------------------------------------------------
now_line(1:5) = [];
%--- parity check ---------------------------------------------------------
now_word = rtcm_parity_check();     % parity bit ����ġ ����ó�� �߰��ؾ���