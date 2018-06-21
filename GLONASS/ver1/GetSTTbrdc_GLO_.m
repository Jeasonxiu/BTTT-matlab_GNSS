function [STT] = GetSTTbrdc_GLO(gs, eph_glo, icol, vec_site, deltat)
%
%function [STT] = GetSTTbrdc_GLO(gs, eph_glo, icol, vec_site, deltat)
%
% DO: Compute the signal transmission time using BRoaDCast ephemeris
%
% <input>   gs: time-epoch given in GPS Week Second
%           eph_glo: EPHemeris array
%           icol: Epoch index for the given <gs>
%           vec_site: Site position
%           deltat

% <output>  STT: signal transmission time [second]
%
% Code by Miso Kim, August 2014

%% ����� �Ű����� ����
CCC = 299792458;
eps = 1.e-10;
MaxIter = 3; % �⺻ 10���� �ӵ� ����� ���� 3������ ������
%% �ݺ� ��� 
stt_0 = 0.075;
for iter = 1:MaxIter
    tau_s = gs - stt_0;
    [sat_pos,sat_vel] = GetSatPosGLO(eph_glo, icol, tau_s, deltat); 
    vec_sat = sat_pos; 
    vec_sat = vec_sat';
    com = norm(vec_sat - vec_site);
    STT = com/CCC;
    if abs(STT - stt_0) < eps
        return
    end
    stt_0 = STT;
end
