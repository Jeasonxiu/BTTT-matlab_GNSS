function [SatPos SatVel] = GetSatVelNC(eph, i_eph, t_GPS)
%
%function [SatVel] = GetSatVelNC(eph, i_eph, t_GPS)
%
% DO: To compute satellite velocity using GetSatPosNC
%
% <input>   eph: ephemeris array
%           i_eph: index of the correspoding epoch close to t_GPS
%           t_GPS: epoch at which position and velocity is required
%   
% <output>  SatVel: Velocity of satellite
%
% Copyright: Kwan-Dong Park, December 9, 2013 @LDEO
%

%% t2 = t1 + dt
dt = 1e-3;
t1 = t_GPS;
t2 = t1 + dt;

pos1 = GetSatPosNC_GC(eph, i_eph, t1);
pos2 = GetSatPosNC_GC(eph, i_eph, t2);
SatPos = pos1;
SatVel = (pos2 - pos1)/ dt;

