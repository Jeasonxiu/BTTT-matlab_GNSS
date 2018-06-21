function [GPGSV] = getGPGSV(filename)
%
% function [GPGSV] = getGPGSV(filename)
%
%   Read the given Logged NMEA file, get GPGSV from NMEA
%   
%   input filename : logged NMEA file
%
%   Example : [] = getGPGSV('NMEA.txt')
%
%   coded by Joonseong Gim, Jan 21, 2016
%   
% GPS Satellites in view
% 
% eg. $GPGSV,3,1,11,03,03,111,00,04,15,270,00,06,01,010,00,13,06,292,00*74
%     $GPGSV,3,2,11,14,25,170,00,16,57,208,39,18,67,296,40,19,40,246,00*74
%     $GPGSV,3,3,11,22,42,067,42,24,14,311,43,27,05,244,00,,,,*4D
% 
% 
%     $GPGSV,1,1,13,02,02,213,,03,-3,000,,11,00,121,,14,13,172,05*67
% 
% 
% 1    = Total number of messages of this type in this cycle
% 2    = Message number
% 3    = Total number of SVs in view
% 4    = SV PRN number
% 5    = Elevation in degrees, 90 maximum
% 6    = Azimuth, degrees from true north, 000 to 359
% 7    = SNR, 00-99 dB (null when not tracking)
% 8-11 = Information about second SV, same as field 4-7
% 12-15= Information about third SV, same as field 4-7
% 16-19= Information about fourth SV, same as field 4-7

fid = fopen(filename,'r');
% fid = fopen('(A)20160111115322.txt');
if fid == -1
    disp('Cannot locate the input file!')
    GPGSV = {};
else
    GPGSV = {};
    fid_out = fopen('GPGSV.txt', 'w');
    GPGSVlist = textscan(fid,'%s');
    
    for j = 1:length(GPGSVlist{1})
        line = cell2mat(GPGSVlist{1}(j));
        if length(line) >= 6
            nmea = line(1:6);
            if nmea == '$GPGSV'
                fprintf(fid_out, '%s \n',line);
                coloum(j,1) = j;
                GPGSV(j,1) = {line};
            end
        end
    end
    if ~isempty(GPGSV)
        GPGSV = GPGSV(find(~cellfun(@isempty,GPGSV)),1);
        fclose(fid_out);    
    end
end    
    





   