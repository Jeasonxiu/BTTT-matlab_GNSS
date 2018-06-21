function [Trop] = Corr_Trop_Saastamoinen(P_rec, P_sat, P)

% Saastamoinen Model(1973)
% ������ ���� ����� �������� ������ ��
% �Էº���: ������ ��ǥ P_rec(m), ���� ��ǥ P_sat(m). ���� P(hPa)
%                ���� ���� ���� ��� 9999�Է�
% ��º���: ����� ��ȣ������(m)

T = 273.16;
e = 0;

% �������� 3���� ��ǥ P_rec�� ����, �浵, Ÿ��ü��� ��ȯ
[lat, lon, h] = xyz2gd(P_rec);

% ���� �Է°��� ������ �ѹݵ��� ���� ���
if P== 9999
    P = (-0.12*h)+1027.6;         % �Էµ� ���� ���� ���� ��� �ѹݵ��� ���� �̿�
end

% ���� E ����
delta = P_sat - P_rec;
[topo] = xyz2topo(delta, lat, lon);
[AzEl] = topo2AzEl(topo);
E = AzEl(2);

Z = 90-E;

Trop = (0.002277/cosd(Z))*(P+(1255/T+0.05)*e-tand(Z)^2);