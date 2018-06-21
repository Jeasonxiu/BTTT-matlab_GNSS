function out = bitcmp2(in, order)
%   OUT = BITCMP2(IN, ORDER)
%   Do : 2�� ����
%  
%   Copyright: taeil Kim, February 20, 2015@INHA University

%% Function
numofbits = @(n_) 2^n_ - 1;     % numofbits(8) = 11111111(2)
%%
if in >= 2^order
    error('������ �ʰ��ϴ� ���� �Է��Ͽ���.')
end

in = double(in);                % ����ȯ uintXX -> double

if in < 0
    sign = 1;
    in = abs(in);               % ����
else
    sign = -1;                  % ���
end

out = sign * bitand( bitcmp(in, order+1) + 1, numofbits(order) );