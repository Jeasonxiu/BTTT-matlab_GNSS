function [out] = IntpSP3(arrSP3, prn, pos, intv) 
%
%function [out] = IntpSP3(arrSP3, prn, pos, intv)
%
% DO: Extract time-tag and interpolated coordinates (pos: 1-X, 2-Y, 3-Z, 4-dt)
%
% <input>   arrSP3: SP3 array
%           prn: PRN number
%           pos: Coordinate indicator
%           intv: Interpolation interval
%
% <out>     out: GPS Week Second and interpolated SP3 <X or Y or Z> 
%               c1(gs), c2(X or Y or Z interpolated)
%
% Copyright: Kwan-Dong Park, November 2013 @LDEO
%
% 주의사항: 출력은 2열로 X/Y/Z 중의 하나만을 보간해서 출력하고 있음

%% 입력변수 테스트
% load 'sp1'; arrSP3 = sp1;
% prn = 32;
% pos = 1; 
%% SP3 행렬에서 몇 번째 칼럼을 사용할 것인가 결정[posArr ~ 3(X), 4(Y), 5(Z)]
posArr = pos + 2;
%% 입력된 arrSP3의 TIME-TAG 추출 및 길이 검토, 길이는 96이어야 함
TT = unique(arrSP3(:,1));
nTT = length(TT);
if nTT ~= 96, error('The number of epochs is not equal to 96'); end
% if nTT ~= 48, error('The number of epochs is not equal to 96'); end
%% 해당 PRN의 좌표 추출: posArr가 3이면 X, 4이면 Y, 5이면 Z, 6이면 dt(시계오차)
id1 = find(arrSP3(:, 2) == prn);
if isempty(id1)
    out = [];
    return
end
    
xPRN = arrSP3(id1, 1);
yPRN = arrSP3(id1, posArr); %: 3=X, 4=Y, 5=Z, 6=dt
%% Lagrange Interpolation의 현재 차수는 9로 세팅되어 있음
nOrder = 9; %: 현재로선 반드시 홀수를 사용해야 함. nRange가 소수점 수가 됨
nPoints = nOrder - 1;
nRange = (nOrder - 1) / 2;
%%
TT0 = TT(1);
TTF = TT0 + 86399; %: TT Final 
nOut = TTF - TT0 + 1;
out = zeros(nOut, 2);
%%
k = 0;
for xint = TT0:intv:TTF %: 보간하려는 시각
    iTT = PickTTrange(TT, xint);
    if iTT - nRange <= 0
        iStart = 1;
        iEnd = nPoints;
    elseif iTT + nRange > nTT
        iStart = nTT - nPoints + 1;
        iEnd = nTT;
    else
        iStart = iTT - nRange + 1;
        iEnd = iTT + nRange;
    end
    %fprintf('%2d: %02d-%02d\n', iTT, iStart, iEnd)
    k = k + 1;    
    xin = xPRN(iStart:iEnd);
    yin = yPRN(iStart:iEnd);
    yint = lagrange(xin, yin, xint);
    out(k,1) = xint;
    out(k,2) = yint;

end
out = out(1:k, :);
%plot(out(:,1), out(:,2), '.b', xPRN, yPRN, 'or')