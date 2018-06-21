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

%% 함수만들기 전 입력변수 테스트
%time = 260200;

%% 입력변수 TT의 길이 계산
nTT = length(TT);

%% 입력된 time이 TT 범위 밖일 때 []을 리턴하도록 함
iTT = [];
if time < TT(1) || time > TT(end)
    disp('** Warning: Input time is out of range [PickTTrange]') 
end

%% time이 포함된 해당 인덱스 찾기
for k = 1:nTT - 1
    if time >= TT(k) && time < TT(k+1);
        iTT = k;
        return;
    else
        iTT = nTT;
    end
end


