function [GPRMC] = getGPRMC(filename)
%
% function [GPRMC] = getGPRMC(filename)
%
%   Read the given Logged NMEA file, get GPRMC from NMEA
%   
%   input filename : logged NMEA file
%
%   Example : [GPRMC] = getGPRMC('NMEA.txt')
%
%   coded by Joonseong Gim, Jan 28, 2016
%   
% RMC - NMEA has its own version of essential gps pvt (position, velocity, time) data. It is called RMC, The Recommended Minimum, which will look similar to:
% 
% $GPRMC,123519,A,4807.038,N,01131.000,E,022.4,084.4,230394,003.1,W*6A
% 
% Where:
%      RMC          Recommended Minimum sentence C
%      123519       Fix taken at 12:35:19 UTC
%      A            Status A=active or V=Void.
%      4807.038,N   Latitude 48 deg 07.038' N
%      01131.000,E  Longitude 11 deg 31.000' E
%      022.4        Speed over the ground in knots
%      084.4        Track angle in degrees True
%      230394       Date - 23rd of March 1994
%      003.1,W      Magnetic Variation
%      *6A          The checksum data, always begins with *


fid = fopen(filename,'r');
% fid = fopen('(A)20160111115322.txt');
% fid = fopen('150707.txt');
% fid = fopen('jprA014a.txt','r');
if fid == -1
    disp('Cannot locate the input file!')
    GPRMC = {};
else
    GPRMC = {};
    fid_out = fopen('GPRMC.txt', 'w');
    GPRMClist = textscan(fid,'%s');

    for j = 1:length(GPRMClist{1})
        line = cell2mat(GPRMClist{1}(j));
        if length(line) >= 6
            nmea = line(1:6);
            if nmea == '$GPRMC' | nmea == '$GNRMC'
                fprintf(fid_out, '%s \n',line);
                coloum(j,1) = j;
                GPRMC(j,1) = {line};
            end
        end
    end
    if ~isempty(GPRMC)
        GPRMC = GPRMC(find(~cellfun(@isempty,GPRMC)),1);
        fclose(fid_out);    
    end
end





   