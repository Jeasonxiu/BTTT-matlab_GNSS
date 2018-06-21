function [tf] = rtcm_preamble(word)
%
%function [tf] = rtcm_preamble(word)
%
%   해당 워드(word)가 한 프레임의 첫번째 워드인지 검사
%
%   Copyright: taeil Kim, February 20, 2015@INHA University

%% Constant
PREAMBLE_MASK= uint32(1069547520);  % 00111111 11000000 00000000 00000000
PREAMBLE     = uint32(427819008);   % 00011001 10000000 00000000 00000000
%%
tf = bitand( word, PREAMBLE_MASK ) == PREAMBLE;