function [yy, mo, dd, hh, mm, ss] = gwgs2date(gw, gs)
% 
% Copyright: Kwan-Dong Park, December 10, 2013 @LDEO
%

%% 함수작성 전 입출력 테스트
% gw = 1770;
% gs = 3540;
%%
jd_DOB = date2jd(1980, 1, 6, 0, 0, 0);             % beginning of GPS week numbering
nWeek = gw;                             
jd = jd_DOB + nWeek*7 + gs/86400;
[yy, mo, dd, hh, mm, ss] = jd2date(jd);
