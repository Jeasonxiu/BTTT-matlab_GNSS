function [dIono] = ionoKlob(al, be, t_GPS, az, el, StaPos)
%
%function [dIono] = ionoKlob(al, be, t_GPS, az, el, StaPos)
%
% do: compute slant ionospheric delay in meters based on Klobuchar model
%
% input: al(1:4): Alphas
%        be(1:4): Betas
%        tGPS
%        azimuth & elevation of the GNSS satellite [degrees] 
%        latitude & longitude of the GNSS site
% output: ionopheric delay in meters
%
StaLLH = xyz2gd(StaPos);            %:conversion of site XYZ to LLH
lat = StaLLH(1);
lon = StaLLH(2);

lat = lat/180; lon = lon/180;       %:conversion to SC(semi-circle) - lat&lon of site
el = el/180; az = az/180;           %:conversion to SC(semi-circle) - azi&ele of sat

CCC = 299792458.;       %:CCC = Speed of Light

Psi = 0.0137/(el + 0.11) - 0.022;   %:Step(1) Earth's central angle between GPS site and IPP
lat_ipp = lat + Psi*cos(az);        %:Step(2) Geodetic latitude of subIPP
if abs(lat_ipp) > 0.416             %:Step(3) IPP conversion
    lat_ipp = 0.416*lat_ipp/abs(lat_ipp);
end

lon_ipp = lon + Psi*sin(az)/cos(lat_ipp);           %:Step(4) Longitude of subIPP
lat_ipp_m = lat_ipp + 0.064*cos(lon_ipp - 1.617);   %:Step(5) Geo-magnetic latitude of subIPP

t = 43200*lon_ipp + t_GPS;                          %:Step(6) t: local time at subIPP; t_GPS: GPS time
if t >= 86400
    t = t - 86400;
elseif t < 0
    t = t + 86400;
end

P = be(1) + be(2)*lat_ipp_m + be(3)*lat_ipp_m^2 + be(4)*lat_ipp_m^3;    %:Step(7) Determination of P
if P < 72000    %:Step(8) Fine-tuning of P
    P = 72000;
end
Q = al(1) + al(2)*lat_ipp_m + al(3)*lat_ipp_m^2 + al(4)*lat_ipp_m^3;    %:Step(9) Determination of Q
if Q < 0        %:Step(10) FIne-tuing of Q
    Q = 0;
end

x = (2*pi*(t - 50400))/P;   %:Step(12) Determination of x
F = 1+16*(0.53 - el)^3;     %:Mapping Function
if abs(x) < 1.57            %:Step(13) Calculation of Ionospheric Delay in meters
    dIono = CCC*F*(5e-9 + Q*(1 - x^2/2 + x^4/4));
else
    dIono = CCC*F*5e-9;
end                         %:Slant TEC = dIono/0.164;
