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
% 1/7/2015: 미소 작성 버전 수정: 특히 회전행렬을 각도3개가 아주 작다는 가정을 해서 연산속도 높임

%% 함수작성 이전 입출력 테스트
% PZ = [1 2 3];
%% 7개 매개변수 c1/2/3(Tx/Ty/Tz [m]), c4/5/6(Rx/Ry/Rz ["]), c7(scale [1+scale])
% p7 = [-0.36 0.08 0.18 0 0 0 0]';  %: REF - Vdovin et al. 2012
p7 = [0.373 -0.186 -0.202 2.3 -3.54 4.21 0.008]';  %: REF - Vdovin et al. 2012
%% 초기화
T = zeros(3,1);
S = zeros(1,1);
R = zeros(3,3);
%% 원점 이동과 축척 정의
T(1:3) = p7(1:3);
S = 1 + p7(7);
%% 회전행렬 - 미소작성 내용을 삭제하고 매우 작은 각도를 가정한 방식으로 교체 
th(1:3) = p7(4:6)*(1/3600)*(pi/180);     %: Arc-Second를 radian으로 변환
R = [ 1       th(3)  -th(2); ...
     -th(3)   1       th(1); ...
      th(2)  -th(1)   1     ];
%% 좌표변환
% WGS = T + S*R*PZ';
% WGS = WGS';
PZ02 = T + S*R*PZ';
PZ02 = PZ02';
%% 김미소 작성사항 삭제
% Rx = [1 0 0; 0 cosd(Ax) sind(Ax); 0 -sind(Ax) cosd(Ax)];
% Ry = [cosd(Ay) 0 -sind(Ay); 0 1 0; sind(Ay) 0 cosd(Ay)];
% Rz = [cosd(Az) sind(Az) 0; -sind(Az) cosd(Az) 0; 0 0 1;];
% Rm = Rz * Ry * Rx;