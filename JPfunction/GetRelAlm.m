function [dRel] = GetRelAlm(alm, prn, gs)
%
% DO: Extract the relativistic effect form Almanac information
%
% <input>   alm: almanac array
%           prn
%           gs: Epoch time to get the relativistic effect
%
% <output>  dRel: Relativity effect
%
% Copyright: Joonseong Gim, January 05, 2017 
%

%% Parameters
e         = alm(find(alm(:,1) == prn),3);
sqtA      = alm(find(alm(:,1) == prn),7);
t_oa      = alm(find(alm(:,1) == prn),4);
M_0       = alm(find(alm(:,1) == prn),10);
%% 상수의 경우 mu만 필요함
mu = 3.986005e14;
%% 계산 시작
A = sqtA^2;
n_0 = sqrt(mu/(A^3));
t_k = gs - t_oa;
M_k = M_0 + n_0*t_k;
E_k = ecce_anom(M_k, e, 3);
%% 상대성 효과 계산
dRel = -4.442807633*10^(-10)*e*sqtA*sin(E_k);
