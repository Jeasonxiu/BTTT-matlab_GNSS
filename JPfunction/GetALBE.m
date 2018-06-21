function [al, be] = GetALBE(eph_file)
%
%function [LeapSec] = GetLeapSec(fileEPH)
%
% Purpose: Read ephemerides file and find AL and BE for Klobuchar model
%
% input: Ephemeris file name
% output: alpha and beta parameters
%
% Copyright: Kwan-Dong Park,13-09-25, LDEO
%

%* Open Epheremides File
fid_eph = fopen(eph_file,'r');

%* Initialization
al = zeros(4,1);
be = zeros(4,1);

%* Counter to set the EOF
ready = 0;

while ~ready
    s = fgetl(fid_eph);
    if length(s) > 71
        if s(61:68) == 'ION ALPH'
            al(1) = str2num(s(4:14));
            al(2) = str2num(s(16:26));
            al(3) = str2num(s(28:38));
            al(4) = str2num(s(40:50));
        end
        if s(61:68) == 'ION BETA'
            be(1) = str2num(s(4:14));
            be(2) = str2num(s(16:26));
            be(3) = str2num(s(28:38));
            be(4) = str2num(s(40:50));
        end
        if s(61:72) =='END OF HEADE'
            ready = 1;
        end
    end
end
