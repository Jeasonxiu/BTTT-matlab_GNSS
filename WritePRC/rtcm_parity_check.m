function [temp] = rtcm_parity_check()
%
%function [temp] = rtcm_parity_check()
%
%   현재 워드(now_word)의 패리티 비트 체크
%
%       pass  : temp = bitxor(temp, DATA_MASK)
%       flase : temp = 0;
%
%   Copyright: taeil Kim, February 20, 2015@INHA University

%% Global
global now_word;
%% Function
oslf_xor  = @(uint_32)  mod(sum(dec2bin(uint_32)-48),2);
%% Constant
PARITY_25 = uint32(3139384448);     % 10 111011 000111 110011 010010 000000
PARITY_26 = uint32(1569692224);     % 01 011101 100011 111001 101001 000000
PARITY_27 = uint32(2932329728);     % 10 101110 110001 111100 110100 000000
PARITY_28 = uint32(1466164864);     % 01 010111 011000 111110 011010 000000
PARITY_29 = uint32(1806824256);     % 01 101011 101100 011111 001101 000000
PARITY_30 = uint32(2340063680);     % 10 001011 011110 101000 100111 000000

P_30_MASK = uint32(1073741824);     % 01 000000 000000 000000 000000 000000
DATA_MASK = uint32(1073741760);     % 00 111111 111111 111111 111111 000000
%%
temp = now_word;
p = uint8(0);

if bitand(temp, P_30_MASK)
    temp = bitxor(temp, DATA_MASK); % d1~d24 위상변환
end

p = bitor( bitshift(p, 1), oslf_xor(bitand(temp, PARITY_25)) );
p = bitor( bitshift(p, 1), oslf_xor(bitand(temp, PARITY_26)) );
p = bitor( bitshift(p, 1), oslf_xor(bitand(temp, PARITY_27)) );
p = bitor( bitshift(p, 1), oslf_xor(bitand(temp, PARITY_28)) );
p = bitor( bitshift(p, 1), oslf_xor(bitand(temp, PARITY_29)) );
p = bitor( bitshift(p, 1), oslf_xor(bitand(temp, PARITY_30)) );

if p ~= bitand( now_word, 63 )      % 63 : 111111
    fprintf('parity bit error\n');
    temp = 0;
end