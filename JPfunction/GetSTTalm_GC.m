function [STT] = GetSTTalm_GC(gs, prn, alm, vec_site)
%
%function [STT] = GetSTTalm_GC(gs, prn, alm, vec_site)
%
% DO: Compute the signal transmission time using BRoaDCast ephemeris
%
% <input>   gs: time-epoch given in GPS Week Second
%           almanac: Almanac array
%           prn
%           vec_site: Site position (1X3)
%
% <output>  STT: signal transmission time [second]
%
% Copyright: Joonseong Gim, January 4, 2017
%
%


%% 상수와 매개변수 정의
CCC = 299792458;
eps = 1.e-10;
MaxIter = 10;

%% 반복 계산 
stt_0 = 0.075;
for iter = 1:MaxIter
    tau_s = gs - stt_0;
    vec_sat = GetSatPos_GC_almanac(alm, prn, tau_s);
    com = norm(vec_sat - vec_site);
    STT = com/CCC;
    if abs(STT - stt_0) < eps
        return
    end
    stt_0 = STT;
end
