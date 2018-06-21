function [SatsEl, highELprn] = PickRSelGLO(gs, LeapSec, Sats, ephglo, Sat_ar, TauC, vec_site)
%
% DO: Find the index of RS(Reference Satellite) with the highest elevation angle
%
% <input>   gs: GPS Week Second
%           Sats: Satellite PRN array
%           eph; Ephemeris array
%           vec_site: Site Vector
%
% <output>  highELprn: Index of RS - PRN with the highest elevation angle
%           SatsEl: List of PRNs with elevation angle
%
% Copyright: Kwan-Dong Park @Jipyong Space, 11/7/2014
%--- Modifications ----
% 11/8/14 �Լ��̸��� RS(��������)�� ���Եǵ��� ������
% 11/23/14 ���� ������ �����ؼ� �����ϵ��� ��. �Ӱ���� ���� ���� ������

%% ������ǥ�� ������ ���� ������ ����� ���� ���浵 ��ȯ
gd = xyz2gd(vec_site); AppLat = gd(1); AppLon = gd(2);
%% �ݺ���� �������� ������ �ִ� ������ �ʱ�ġ
highEL = 0.;
%% �������� ���� �� ��� ���� �ʱ�ȭ
NoSats = length(Sats);
SatsEl = NaN(NoSats, 3);
SatsEl(:, 1) = gs;
%% �ݺ���� �������� �ִ���� ���� ����
for kS = 1:NoSats;
    prn = Sats(kS);
    
    tc = gs - LeapSec;
    icol=PickEPH_GLO2(ephglo, prn, tc);
    
    TauN=ephglo(icol,12); GammaN=ephglo(icol,13); %: tau & gamma �ð���� ������ ���
    ch_num=ephglo(icol,16); %: channel number ������ ������ ���
    
    % ��ȣ���޽ð� ���
    STT = GetSTTbrdcGLO2(Sat_ar,gs,prn,vec_site');
    % LeapSecond & ��ȣ���� �ð��� ������ ���� ��ġ ����
    [vec_sat, SatVel] = SatPosLS_STT(Sat_ar,gs,prn,LeapSec,STT,TauC);
    vec_rho = vec_sat - vec_site;
    
    [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
    
    
    if el > highEL
        highEL = el;
        highELprn = kS;
    end
    
    SatsEl(kS, 2) = prn;
    SatsEl(kS, 3) = el;
    
    %     if el < 15
    %         disp([gs prn el])
    %     end
    %     fprintf(1, '%8d %3d %6.2f \n', gs, prn, el)
end