function [out] = IntpBRDC_glo1e(eph_glo, prn, tt)
%
%function [out] = IntpBRDC_glo1e(eph_glo, prn, tt)
%
% DO: Extract time-tag and interpolated coordinates
%
% <input>   eph_glo: glonass eph array using ReadEPH_glo_kdpark function
%           prn: PRN number
%           tt: Time-Tag to interpolate
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

if nTT ~= 48, error('The number of epochs is not equal to 48'); end
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

%% Interpolation 시작
tint = tt;

%% 입력 tt 가 eph_glo Range 안에 존재할 경우에만 동작
%% eph_glo 첫 Time-Tag 이전, 마지막 Time-Tag 보다 30분 이후 경우 동작 안함
if tint < TT(1) || tint > TT(end)
    disp('out of range');
    out = [];
    return
else
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
    
    tin = PRNt(iStart:iEnd);
    xin = PRNx(iStart:iEnd);
    yin = PRNy(iStart:iEnd);
    zin = PRNz(iStart:iEnd);
    xint = lagrange(tin, xin, tint);
    yint = lagrange(tin, yin, tint);
    zint = lagrange(tin, zin, tint);
    out(1,1) = tint;
    out(1,2) = xint;
    out(1,3) = yint;
    out(1,4) = zint;
end

