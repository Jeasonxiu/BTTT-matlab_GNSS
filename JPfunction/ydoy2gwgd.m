function [gw, gd] = ydoy2gwgd(yy, doy)
%
%function [gw, gd] = ydoy2gwgd(yy, doy)
%
% DO: Convert Year & DOY to GPS Week and GPS Week Day Number
%
% <input>   yy:  Year          eg) 14
%           doy: Day of Year   eg) 100
%
% <output>  gw/gd: GPS Week Number and GPS Week Day Number
%
% Copyright: Kwan-Dong Park, October 10, 2014@Jipyong Space
%            - 김미소 버전(10/9/14)을 약간 수정함
%

%% 2자리 년도를 4자리 년도로 변경
if yy < 80
    y4 = 2000 + yy;
else
    y4 = 1900 + yy;
end
%% 윤년 행렬 설정 - 10년후에 변경 필요
switch (y4)
    case {1980,1984,1988,1992,1996,2000,2004,2008,2012,2016,2020,2024}
        days = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    otherwise
        days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
end
%% Calendar Month/Day 결정
tdays = 0;
for i = 1: 12
    tdays = tdays + days(i);
    if doy - tdays < 0
        break;
    end
end
mo = i;
dd = doy - sum(days(1:i-1));
%% JD 계산: date2jd 응용 h/m/s를 12/0/0으로 설정하고 계산
JD = 367*y4 - floor(7*(y4 + floor((mo + 9)/12))/4) + floor(275*mo/9) + dd + 1721013.5 + 0.5;
%% JD 기반으로 GPS Week과 GPS Week Day 계산: jd2gwgs 활용
[gw, gs] = jd2gwgs(JD);
gd = floor(gs/86400);