function out = bitcmp2(in, order)
%   OUT = BITCMP2(IN, ORDER)
%   Do : 2의 보수
%  
%   Copyright: taeil Kim, February 20, 2015@INHA University

%% Function
numofbits = @(n_) 2^n_ - 1;     % numofbits(8) = 11111111(2)
%%
if in >= 2^order
    error('범위를 초과하는 수를 입력하였음.')
end

in = double(in);                % 형변환 uintXX -> double

if in < 0
    sign = 1;
    in = abs(in);               % 음수
else
    sign = -1;                  % 양수
end

out = sign * bitand( bitcmp(in, order+1) + 1, numofbits(order) );