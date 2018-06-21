function [STT] = GetSTTbrdc(gs, prn, eph, vec_site)
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
% Copyright: Kwan-Dong Park, January 17, 2014
%
% RotSatPos�� �ʿ��� ������ ��ȣ���޽ð��� �ݺ������� �����
%
%-- Modifications --
% 8/10/14: �Ű����� �Է¼����� GetSTTsp3�� �����ϰ� ����; (gs, prn, eph, vec_site)

%% ����� �Ű����� ����
CCC = 299792458;
eps = 1.e-10;
MaxIter = 10;
%% icol���� - 8/10/15 ��������
icol = PickEPH(eph, prn, gs);
%% �ݺ� ��� 
stt_0 = 0.075;
for iter = 1:MaxIter
    tau_s = gs - stt_0;
    vec_sat = GetSatPosNC(eph, icol, tau_s);
    com = norm(vec_sat - vec_site);
    STT = com/CCC;
    if abs(STT - stt_0) < eps
        return
    end
    stt_0 = STT;
end
