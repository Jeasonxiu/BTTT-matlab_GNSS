function [out] = IntpBRDC_glo(eph_glo, prn, intv, LeapSec) 
%
%function [out] = IntpBRDC_glo(eph_glo, prn, intv, LeapSec) 
%
% DO: Extract time-tag and interpolated coordinates 
%
% <input>   eph_glo: glonass eph array using ReadEPH_glo_kdpark function
%           prn: PRN number
%           intv: Interpolation interval(second)
%           LeapSec
%
% <out>     out: GPS Week Second and interpolated BRDC [gs, X, Y, Z]
%
% Copyright: Kwan-Dong Park, November 2013 @LDEO
%

%% 입력변수 테스트
% load 'sp1'; arrSP3 = sp1;
% prn = 32;
% pos = 1; 

%% 입력된 eph_glo의 TIME-TAG 추출 및 길이 검토, 길이는 48이어야 함
TT = unique(eph_glo(:,1));
nTT = length(TT);

if nTT ~= 48, error('The number of epochs is not equal to 96'); end
%% 해당 PRN의 좌표 추출
id1 = find(eph_glo(:, 2) == prn);
if isempty(id1)
    out = [];
    return
end

PRNt = eph_glo(id1, 1); % gps second
PRNx = eph_glo(id1, 3); % x coordinate
PRNy = eph_glo(id1, 4); % y coordinate
PRNz = eph_glo(id1, 5); % z coordinate
%% Lagrange Interpolation의 현재 차수는 9로 세팅되어 있음
nOrder = 9; %: 현재로선 반드시 홀수를 사용해야 함. nRange가 소수점 수가 됨
nPoints = nOrder - 1;
nRange = (nOrder - 1) / 2;
%%
TT0 = TT(1);
TTF = TT0 + 84899; %: TT Final 
nOut = TTF - TT0 + 1;
out = zeros(nOut, 2);
%% Interpolation 시작
k = 0;
for tint = TT0:intv:TTF %: 보간하려는 시각
    iTT = PickTTrange(TT, tint);
    if iTT - nRange <=0
        iStart = 1;
        iEnd = nPoints;
    elseif iTT + nRange > nTT
        iStart = nTT - nPoints + 1;
        iEnd = nTT;
    else
        iStart = iTT - nRange + 1;
        iEnd = iTT + nRange;
    end
    
    k = k + 1;    
    tin = PRNt(iStart:iEnd);
    xin = PRNx(iStart:iEnd);
    yin = PRNy(iStart:iEnd);
    zin = PRNz(iStart:iEnd);
    xint = lagrange(tin, xin, tint-LeapSec);
    yint = lagrange(tin, yin, tint-LeapSec);
    zint = lagrange(tin, zin, tint-LeapSec);
    out(k,1) = tint;
    out(k,2) = xint;
    out(k,3) = yint;
    out(k,4) = zint;

end
out = out(1:k, :);
