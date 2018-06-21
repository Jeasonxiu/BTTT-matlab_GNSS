function [STT] = GetSTTsp3(gs, prn, sp3, vec_site)
%
%function [STT] = GetSTTsp3(gs, prn, sp3, vec_site)
%
% DO: Compute the signal transmission time using SP3 ephemeris
%
% <input>   gs: time-epoch given in GPS Week Second
%           prn: PRN
%           sp3: SP3 ephemeris array
%           vec_site: Site position
%
% <output>  STT: signal transmission time [second]
%
% Copyright: Kwan-Dong Park, January 21, 2014
%   - Slightly modified from GetSTTbrdc
%
% RotSatPos에 필요한 것으로 신호전달시간을 반복적으로 계산함
%

%% 상수와 매개변수 정의
CCC = 299792458;
eps = 1.e-10;
MaxIter = 3;
%% 반복 계산 
stt_0 = 0.075;
for iter = 1:MaxIter
    tau_s = gs - stt_0;
    out4 = IntpSP3e1(sp3, prn, tau_s); 
    vec_sat = out4(2:4)';
    com = norm(vec_sat - vec_site);
    STT = com/CCC;
    if abs(STT - stt_0) < eps
        return
    end
    stt_0 = STT;
end
