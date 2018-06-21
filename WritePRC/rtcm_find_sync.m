function [] = rtcm_find_sync()
%
%function [] = rtcm_find_sync()
%
%   아무 프레임의 첫번째 워드를 찾음
%
%   Copyright: taeil Kim, February 20, 2015@INHA University

%% Global
global now_line;
global next_line;
global now_word;
global sync;
%%
%--- now_word << 6 bit ----------------------------------------------------
if length(now_line) < 1
    now_line = [now_line next_line(1)];
    next_line(1) = [];
end
u_char = uint32( now_line(1) );
now_line(1) = [];           % 사용한 byte 삭제
for i=1:6
    now_word = bitor( bitshift(now_word, 1), bitget(u_char, i) );
end                         % reverse
%--- 프레임 시작 찾기 ------------------------------------------------------
check = rtcm_parity_check();
if check
    if rtcm_preamble(check)
        sync = 1;
        now_word = check;
    end
end