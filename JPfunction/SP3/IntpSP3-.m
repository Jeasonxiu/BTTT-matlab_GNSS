function [out] = IntpSP3(arrSP3, prn, pos) 
%
% DO: Extract time-tag and interpolated coordinates (pos: 1-X, 2-Y, 3-Z, 4-dt)
%

%% 입력변수 테스트
% load 'sp1'; arrSP3 = sp1;
% prn = 32;
% pos = 1; 
%%
posArr = pos + 2;
%% 입력된 arrSP3의 TIME-TAG 추출 및 길이 검토, 길이는 96이어야 함
TT = unique(arrSP3(:,1));
nTT = length(TT);
if nTT ~= 96, error('The number of epochs is not equal to 96'); end  % 96 = 24 * 4 (15min interval)
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
% TTF = TT0 + (86400-1); %: TT Final
TTF = TT0 + (86400-(15*60)); %: TT Final
nOut = TTF - TT0 + 1;
out = zeros(nOut, 2);
%%
k = 1;
for xint = TT0:1:TTF %: 보간하려는 시각
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
    
    xin = xPRN(iStart:iEnd);
    yin = yPRN(iStart:iEnd);
    yint = lagrange(xin, yin, xint);
    out(k,1) = xint;
    out(k,2) = yint;
    k = k + 1;
end
%plot(out(:,1), out(:,2), '.b', xPRN, yPRN, 'or')