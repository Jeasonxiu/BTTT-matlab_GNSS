function [PZ02] = PZ112PZ02(PZ)
%
%function [WGS] = PZ2WGS(PZ)
% 
% DO: Coordinate Transformation from PZ90.02 to WGS84.
%
% <input>   PZ:  XYZ in PZ90  [1 X 3]
%
% <output>  WGS: XYZ in WGS84 [1 X 3]
%
% Copyright: Miso Kim & Kwan-Dong Park, January 7th, 2015
%
%--- Modifications --
% 1/7/2015: �̼� �ۼ� ���� ����: Ư�� ȸ������� ����3���� ���� �۴ٴ� ������ �ؼ� ����ӵ� ����

%% �Լ��ۼ� ���� ����� �׽�Ʈ
% PZ = [1 2 3];
%% 7�� �Ű����� c1/2/3(Tx/Ty/Tz [m]), c4/5/6(Rx/Ry/Rz ["]), c7(scale [1+scale])
% p7 = [-0.36 0.08 0.18 0 0 0 0]';  %: REF - Vdovin et al. 2012
p7 = [0.373 -0.186 -0.202 2.3 -3.54 4.21 0.008]';  %: REF - Vdovin et al. 2012
%% �ʱ�ȭ
T = zeros(3,1);
S = zeros(1,1);
R = zeros(3,3);
%% ���� �̵��� ��ô ����
T(1:3) = p7(1:3);
S = 1 + p7(7);
%% ȸ����� - �̼��ۼ� ������ �����ϰ� �ſ� ���� ������ ������ ������� ��ü 
th(1:3) = p7(4:6)*(1/3600)*(pi/180);     %: Arc-Second�� radian���� ��ȯ
R = [ 1       th(3)  -th(2); ...
     -th(3)   1       th(1); ...
      th(2)  -th(1)   1     ];
%% ��ǥ��ȯ
% WGS = T + S*R*PZ';
% WGS = WGS';
PZ02 = T + S*R*PZ';
PZ02 = PZ02';
%% ��̼� �ۼ����� ����
% Rx = [1 0 0; 0 cosd(Ax) sind(Ax); 0 -sind(Ax) cosd(Ax)];
% Ry = [cosd(Ay) 0 -sind(Ay); 0 1 0; sind(Ay) 0 cosd(Ay)];
% Rz = [cosd(Az) sind(Az) 0; -sind(Az) cosd(Az) 0; 0 0 1;];
% Rm = Rz * Ry * Rx;