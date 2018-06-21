function [iTT] = PickTTrange(TT, time)

%
%function [icol] = PickTTrange(TT, time)
%
% DO: Find the index of the range where the input 'time' resides
%
% <input>   TT: Array of time tags
%           time: tk
%
% <output>  iTT: index of TT where 'time' is located
%
% Copyright: Kwan-Dong Park@LDEO, October 26, 2013
%

%% �Լ������ �� �Էº��� �׽�Ʈ
%time = 260200;

%% �Էº��� TT�� ���� ���
nTT = length(TT);

%% �Էµ� time�� TT ���� ���� �� []�� �����ϵ��� ��
iTT = [];
if time < TT(1) || time > TT(end)
    disp('** Warning: Input time is out of range [PickTTrange]') 
end

%% time�� ���Ե� �ش� �ε��� ã��
for k = 1:nTT - 1
    if time >= TT(k) && time < TT(k+1);
        iTT = k;
        return;
    else
        iTT = nTT;
    end
end


