function [VelDop] = GetVelDop(gs, prns, dops, eph, vec_site)
%
%function [VelDop] = GetVelDop(gs, prns, dops, eph, vec_site)
%
% DO: Compute velocity based on Doppler measurements
%
% <input>   gs: GPS week second
%           prns: array of prns observed at the epoch
%           dops: doppler measurements at the epoch correspong to the prn
%           eph: ephemeris array read using ReadEPH
%           vec_site: receiver position at the epoch
%
% <output>  VelDop: Velocity computed from Doppler observables
%
% Copyright: Kwan-Dong Park@Jipyong Space, 7/12/2015
%
%% �ʿ��� ��� �� �Ű����� ����
CCC = 299792458.; % CCC = Speed of Light [m/s]
f1 = 1575.42e6;   % L1 ���ļ�
lmbd1 = CCC/f1;   % L1 ����
%% �ʱⰪ 4���� ��� 0���� �����ص� ū ���� ���� ������ �ľ��Ͽ���
x = [0 0 0 0]; x = x';
vel_site = x(1:3)';
%% ���Թ����Ŀ� �ʿ��� ��� �ʱ�ȭ
HTH = zeros(4,4);
HTy = zeros(4,1);
%% ���� �ּ������� ���� ����
NoSats = length(prns);
for kS = 1:NoSats
    prn = prns(kS);
    obs = -dops(kS);
    
    icol = PickEPH(eph, prn, gs);
    [vec_sat,vel_sat] = GetSatVelNC(eph, icol, gs);
    
    vec_site2sat = vec_sat - vec_site;
    vel_site2sat = vel_sat - vel_site;
    
    com = dot(vec_site2sat/norm(vec_site2sat), vel_site2sat) + x(4);
    y = lmbd1*obs - com;
    
    H = [ -vec_site2sat/norm(vec_site2sat) 1];
    
    HTH = HTH + H'*H;
    HTy = HTy + H'*y;
end
%% ���������̹Ƿ� �ݺ� ��� �ʿ� ����
x = x + inv(HTH)*HTy;
VelDop = x(1:3)';
