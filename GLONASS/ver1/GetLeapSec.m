function [LeapSec] = GetLeapSec(fid_eph)
%
%function [LeapSec] = GetLeapSec(fid_eph)
%
%   Read ephemerides file and find the leap seconds
%       - Usually called by ReadEPH
%
%   Originally coded by an Unknown person,
%   Modified by Kwan-Dong Park, July 23, 2006, Kookmin University, SSOA

fid_eph = fopen(fid_eph,'r');
%* Set LeapSec = 0 in case the header file doesn't have the leap seconds
LeapSec = 0;

%* Counter to set the EOF
ready = 0;

while ~ready
    s = fgetl(fid_eph);
    if length(s) > 71
        if s(61:72) == 'LEAP SECONDS'
            LeapSec = str2num(s(1:6));
        end
        if s(61:72)=='END OF HEADE'
            ready = 1;
        end
    end
end
