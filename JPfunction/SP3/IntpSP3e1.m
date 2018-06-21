function [out] = IntpSP3e1(arrSP3, prn, tt) 
%
%function [out] = IntpSP3e1(arrSP3, prn, tt)
%
% DO: Extract time-tag and interpolated coordinates X, Y, and Z for 
%     just one epoch (referred to as 'e1')
%
% <input>   arrSP3: SP3 array
%           prn: PRN number
%           tt: Time-Tag to interpolate
%
% <out>     out: Interpolated X, Y, Z
%               c1(X), c2(Y), c3(Z)
%
% Copyright: Kwan-Dong Park, December 2013 @LDEO
%
% 주의: IntpSP3e1은 한 epoch에 대해서, IntpSP3는 단번에 모든 epoch에 대해 보간
% -Modifications-
% 1/20/2014: tt가 주어진 시간범위(0초) 이전일 때 문제가 생겨서, 조건문을 추가함 tint<TT(1)
%           - 그러나, 이 방법은 신호전달시간 계산과 같은 근접한 시각에만 사용하는 것이 좋을 것임
% 11/30/2014: 출력out을 row 벡터로 변경
% 8/8/2015: 출력 out에서 시간 삭제 - 타 함수와 같이 위성궤도만 출력! 기존 버전 _old

%% 입력변수 테스트
% load 'SP3'; arrSP3 = sp3;
% prn = 1;
% pos = 1;
% tt =  55;
%% 입력된 arrSP3의 TIME-TAG 추출 및 길이 검토, 길이는 96이어야 함
TT = unique(arrSP3(:,1));
nTT = length(TT);
if nTT ~= 96, error('The number of epochs is not equal to 96'); end
%% 해당 PRN의 좌표 추출: posArr가 3이면 X, 4이면 Y, 5이면 Z, 6이면 dt(시계오차)
id1 = find(arrSP3(:, 2) == prn);
if isempty(id1)
    out = [];
    return
end
    
tPRN = arrSP3(id1, 1);
xPRN = arrSP3(id1, 3); %: 3=X, 4=Y, 5=Z, 6=dt
yPRN = arrSP3(id1, 4);
zPRN = arrSP3(id1, 5); 
%% Lagrange Interpolation의 현재 차수는 9로 세팅되어 있음
nOrder = 9; %: 현재로선 반드시 홀수를 사용해야 함. nRange가 소수점 수가 됨
nPoints = nOrder - 1;
nRange = (nOrder - 1) / 2;
%%
out = zeros(1, 3);  %: 출력행렬 크기 사전 설정
tint = tt;          %: 보간하려는 시각
iTT = PickTTrange(TT, tint);

if tint < TT(1) || iTT - nRange <= 0 %: 시계열 이전 시각이 입력될 경우에 대비
    iStart = 1;
    iEnd = nPoints;
elseif iTT + nRange > nTT
    iStart = nTT - nPoints + 1;
    iEnd = nTT;
else
    iStart = iTT - nRange + 1;
    iEnd = iTT + nRange;
end
%%
tin = tPRN(iStart:iEnd);
xin = xPRN(iStart:iEnd);
yin = yPRN(iStart:iEnd);
zin = zPRN(iStart:iEnd);
xout = lagrange(tin, xin, tint);
yout = lagrange(tin, yin, tint);
zout = lagrange(tin, zin, tint);
out(1) = xout;
out(2) = yout;
out(3) = zout;
