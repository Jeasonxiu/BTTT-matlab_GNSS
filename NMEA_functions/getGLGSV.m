function [GLGSV] = getGLGSV(filename)
%
% function [GLGSV] = getGLGSV(filename)
%
%   Read the given Logged NMEA file, get GLGSV from NMEA
%   
%   input filename : logged NMEA file
%
%   Example : [] = getGLGSV('NMEA.txt')
%
%   coded by Joonseong Gim, June 1, 2016
%   
% GPS Satellites in view
% 
% eg. $GLGSV,3,1,10,66,00,015,,67,00,062,36,73,89,233,48,74,35,329,50*61
%     $GLGSV,3,2,10,80,30,151,45,82,14,039,50,83,65,005,52,84,47,244,48*61
%     $GLGSV,3,3,10,85,00,229,,,,,51*54
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
    GLGSV = {};
else
    GLGSV = {};
    fid_out = fopen('GLGSV.txt', 'w');
    GLGSVlist = textscan(fid,'%s');
    
    for j = 1:length(GLGSVlist{1})
        line = cell2mat(GLGSVlist{1}(j));
        if length(line) >= 6
            nmea = line(1:6);
            if nmea == '$GLGSV'
                fprintf(fid_out, '%s \n',line);
                coloum(j,1) = j;
                GLGSV(j,1) = {line};
            end
        end
    end
    if ~isempty(GLGSV)
        GLGSV = GLGSV(find(~cellfun(@isempty,GLGSV)),1);
        fclose(fid_out);    
    end
end    
    





   