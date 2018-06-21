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
% 회전, 스케일 무시해도 mm 수준의 차이임. 주석처리 참고
% Copyright: Jihye Won, December 26th, 2016
% REF - CLONASS ICD 5.1 & PZ-90.11 reference document(2014)

%% 옵셋만 고려하는 버전
WGS =  PZ' - [-0.013 0.106 0.022]'; % PZ90.11 보정값
WGS = WGS';

%% 회전, 스케일까지 고려하는 버전 (PZ-90.11 reference document,2014)

% WGS = (1 + 0.008*10^-6)*[1 -0.65067*10^-6 0.01716*10^-6; 
%                        0.65067*10^-6 1 0.01115*10^-6;
%                       -0.01716*10^-6 -0.01115*10^-6 1] * PZ' + [0.013; -0.106; -0.022];
%WGS=WGS';
end