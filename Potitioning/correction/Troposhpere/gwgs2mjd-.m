function [mjd] = gwgs2mjd(gw, gs)
% 
% Copyright: Kwan-Dong Park, December 15, 2013 @LDEO
%

%% 함수작성 전 입출력 테스트
% gw = 1770;
% gs = 3540;
%%
jd_DOB = date2jd(1980, 1, 6, 0, 0, 0);  %: Date of Birth: 01-06-1980                        
jd = jd_DOB + gw*7 + gs/86400;
mjd = jd - 2400000.5;
