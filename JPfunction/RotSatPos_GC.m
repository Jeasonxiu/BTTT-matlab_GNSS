function [SatPos] = RotSatPos_GC(SatPos, STT, prn)
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
if prn < 200
    ome_e = 7.2921151467e-5; %: Earth Rotation Rate [rad/s]
else
    ome_e = 7.2921150e-5; %: Earth Rotation Rate [rad/s]
end
rota = ome_e * STT;
R_e = [  cos(rota) sin(rota) 0;
        -sin(rota) cos(rota) 0;
            0         0      1];
SatPos = R_e * SatPos';
