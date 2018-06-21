function [Trop] = Corr_Trop_Hopfield(P_rec, P_sat, P)

% Hopfield Model(1969)
% ��ǥ���� ���������� �̿��Ͽ� ���� ���� ������ ��ȭ�� ��Ÿ�� ��
% �Էº���: ������ ��ǥ P_rec(m), ���� ��ǥ P_sat(m). ���� P(hPa)
%                ���� ���� ���� ��� 9999�Է�
% ��º���: ����� ��ȣ������

T = 273.16;
e = 0;

% �������� 3���� ��ǥ P_rec�� ����, �浵, Ÿ��ü��� ��ȯ
[lat, lon, h] = xyz2gd(P_rec);

% ���� �Է°��� ������ �ѹݵ��� ���� ���
if P == 9999
    P = (-0.12*h)+1027.6;        % �Էµ� ������ ���� ��� �ѹݵ��� ���� �̿�
end

% ���� E ����
delta = P_sat - P_rec;
[topo] = xyz2topo(delta, lat, lon);
[AzEl] = topo2AzEl(topo);
E = AzEl(2);

Trop_d = (10^(-6)/5)*(77.64/(sind(sqrt(E^2+6.25))))*P/T*(40136+148.72*(T-273.16));
Trop_w = (10^(-6)/5)*((-12.96*T + 3.718*10^5)/(sind(sqrt(E^2+2.25)))) * (e/(T^2)) * 11000;
Trop = Trop_d + Trop_w;