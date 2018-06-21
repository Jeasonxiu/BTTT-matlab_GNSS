function [WGS] = PZ2WGS(PZ)
%
%function [WGS] = PZ2WGS(PZ)
% 
% DO: Coordinate Transformation from PZ90.11 to WGS84.
%
% <input>   PZ:  XYZ in PZ90  [1 X 3]
%
% <output>  WGS: XYZ in WGS84 [1 X 3]
%
% Copyright: Joonseong GIM, August 25th, 2016
%

%% �Լ��ۼ� ���� ����� �׽�Ʈ
% PZ = [1 2 3];
%% 7�� �Ű����� c1/2/3(Tx/Ty/Tz [m]), c4/5/6(Rx/Ry/Rz ["]), c7(scale [1+scale])
% p7 = [-0.36 0.08 0.18 0 0 0 0]';  %: REF - Vdovin et al. 2012
p7 = [0.013 -0.106 -0.022 0.0023 -0.00354 0.00421 0.000000008]';  %: REF - PZ-90.11 2014 Reference Document
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
WGS = T + S*R*PZ';
WGS = WGS';
