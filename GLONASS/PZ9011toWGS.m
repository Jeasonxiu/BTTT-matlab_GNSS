function [WGS] = PZ9011toWGS(PZ)
% 
% function [newx] = PZ9011toWGS(x)
%  
% DO: Coordinate Transformation from PZ90.11 to WGS84(G1150).
%
% <input>   PZ:  XYZ in PZ90.11  [1 X 3]
%
% <output>  WGS: XYZ in WGS84 [1 X 3]
%
% Copyright: Jihye Won & Hyun-Woo Tae & Joonseong GIM, December 26th, 2016
% ȸ��, ������ �����ص� mm ������ ������. �ּ�ó�� ����
% Copyright: Jihye Won, December 26th, 2016
% REF - CLONASS ICD 5.1 & PZ-90.11 reference document(2014)

%% �ɼ¸� ����ϴ� ����
WGS =  PZ' - [-0.013 0.106 0.022]'; % PZ90.11 ������
WGS = WGS';

%% ȸ��, �����ϱ��� ����ϴ� ���� (PZ-90.11 reference document,2014)

% WGS = (1 + 0.008*10^-6)*[1 -0.65067*10^-6 0.01716*10^-6; 
%                        0.65067*10^-6 1 0.01115*10^-6;
%                       -0.01716*10^-6 -0.01115*10^-6 1] * PZ' + [0.013; -0.106; -0.022];
%WGS=WGS';
end