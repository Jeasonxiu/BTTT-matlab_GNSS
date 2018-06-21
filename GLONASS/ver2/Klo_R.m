function [delta_I_freq] = Klo_R(app_xyz, Alpha, Beta, t_GPS, sat_pos,freq_num)
% klobuchar model for GLONASS
% function [delta_I_freq] = Klo_R(app_xyz, Alpha, Beta, t_GPS, sat_pos,freq_num)
%
% <input>   app_xyz:
%           Alpha:
%           Beta: 
%           t_GPS:
%           sat_pos:
%           freq_num: 
%
% <output>  delta_I_freq: Ionospheric delay [m]
%
% Copyright: Hye-In Kim, October 10, 2014

app_llh = xyz2gd(app_xyz);
lat = app_llh(1);
lon = app_llh(2);
delta_xyz = sat_pos-app_xyz;
delta_nev = xyz2topo(delta_xyz, lat, lon);
AzEl = topo2AzEl(delta_nev);
Az = AzEl(1);
El = AzEl(2);

lat = lat/180;
lon = lon/180;

El = El/180;
Az = Az/180;

c = 299792458.;

% 1. calculted the earth-centred angle
Psi = (0.0137/(El+0.11))-0.022;

% 2. compute the latitude of Ionospheric Pierce Point(IPP)
lat_ipp = lat+(Psi*cos(Az));
if abs(lat_ipp) > 0.416
    lat_ipp = (0.416*lat_ipp)/abs(lat_ipp);
end

% 3. compute the longitude of the IPP
lon_ipp = lon+(Psi*(sin(Az)/cos(lat_ipp)));

% 4. find the geomagnetic latitude of the IPP
lat_ipp_m = lat_ipp+(0.064*cos(lon_ipp-1.167));

% 5. find the local time at the IPP
t = (43200*lon_ipp)+t_GPS;
if t >= 86400
    t = t-86400;
elseif t < 0
    t = t+86400;
end

% 6. compute the amplitude of ionospheric delay
Q = Alpha(1)+(Alpha(2)*lat_ipp_m)+(Alpha(3)*(lat_ipp_m)^2)+(Alpha(4)*(lat_ipp_m)^3);
if Q < 0
    Q = 0;
end

% 7. compute the period of ionospheric delay
P = Beta(1) + (Beta(2)*lat_ipp_m) + (Beta(3)*(lat_ipp_m)^2) + (Beta(4)*(lat_ipp_m)^3);
if P < 72000
    P = 72000;
end

% 8. compute the phase of ionospheric delay
x = (2*pi*(t-50400))/P;

% 9. compute the slant factor 
F = 1+16*(0.53-El)^3;

% 10. compute the ionospheric time delay
if abs(x) < 1.57 
    delta_I = c*F*((5*10^-9)+Q*(1-(x^2/2)+(x^4/4)));
else
    delta_I = c*F*(5*10^-9);
end

freq_GL1=1575.42;   % GPS L1 freq
freq_target=1602 + freq_num*.5625;
delta_I_freq=(freq_GL1/freq_target)^2 * delta_I;

TECU=delta_I/0.164;
