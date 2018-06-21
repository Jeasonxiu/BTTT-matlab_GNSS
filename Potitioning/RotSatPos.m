function [SatPos] = RotSatPos(SatPos, STT)
%
%function [SatPos] = RotSatPos(SatPos, STT)
%
% DO: Rotate the satellite position during the signal transmission time
%
% <input>   SatPos: Satellite position 1x3 [m]
%           STT: Signal transmission time
%    
% <output>  SatPos: Rotated satellite position [m]
%
% REF: Misra and Enge, p202
%
% Copyright: Kwan-Dong Park, January 17, 2014
%
%-- Modifications --
% 8/10/14 출력되는 위성좌표를 1X3 열벡터(row vector)로 변환
%
ome_e = 7.2921151467e-5; %: Earth Rotation Rate [rad/s]
rota = ome_e * STT;
R_e = [  cos(rota) sin(rota) 0;
        -sin(rota) cos(rota) 0;
            0         0      1];
SatPos = R_e * SatPos';
%% 출력되는 위성좌표를 1X3 열벡터(row vector)로 변환
SatPos = SatPos';
