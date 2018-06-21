function [STT] = GetSTTbrdc_GLO_ver3(gs, SatPosArr_before, jcol, vec_site, deltat)
%
%function [STT] = GetSTTbrdc(gs, eph, icol, vec_site)
%
% DO: Compute the signal transmission time using BRoaDCast ephemeris
%
% <input>   gs: time-epoch given in GPS Week Second
%           eph: EPHemeris array
%           icol: Epoch index for the given <gs>
%           vec_site: Site position
%
% <output>  STT: signal transmission time [second]
%
% RotSatPos에 필요한 것으로 신호전달시간을 반복적으로 계산함
%

%% 상수와 매개변수 정의
CCC = 299792458;
eps = 1.e-10;
MaxIter = 10;
%% 반복 계산 
stt_0 = 0.075;
for iter = 1:MaxIter
    tau_s = gs - stt_0;
    vec_sat = GetSatPosGLO_ver3(SatPosArr_before, jcol, tau_s, deltat); vec_sat = vec_sat';
    com = norm(vec_sat - vec_site);
    STT = com/CCC;
    if abs(STT - stt_0) < eps
        return
    end
    stt_0 = STT;
end
