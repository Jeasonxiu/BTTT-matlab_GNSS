function [GPGSA] = getGPGSA(filename)
%
% function [GPGSA] = getGPGSA(filename)
%
%   Read the given Logged NMEA file, get GPGSV from NMEA
%   
%   input filename : logged NMEA file
%
%   Example : [] = getGPGSA('NMEA.txt')
%
%   coded by Joonseong Gim, Jan 22, 2016
%   
% GPS DOP and active satellites
% 
% eg1. $GPGSA,A,3,,,,,,16,18,,22,24,,,3.6,2.1,2.2*3C
% eg2. $GPGSA,A,3,19,28,14,18,27,22,31,39,,,,,1.7,1.0,1.3*35
% 
% 



fid = fopen(filename,'r');
% fid = fopen('ubx2_150707_10m_nmea.txt');

if fid == -1
    disp('Cannot locate the input file!')
    GPGSA = {};
else
    GPGSA = {};
    fid_out = fopen('GPGSV.txt', 'w');
    GPGSAlist = textscan(fid,'%s');
    
    for j = 1:length(GPGSAlist{1})
        line = cell2mat(GPGSAlist{1}(j));
        if length(line) >= 6
            nmea = line(1:6);
            if nmea == '$GPGSA'
                fprintf(fid_out, '%s \n',line);
                coloum(j,1) = j;
                GPGSA(j,1) = {line};
            end
        end
    end
    if ~isempty(GPGSA)
        GPGSA = GPGSA(find(~cellfun(@isempty,GPGSA)),1);
        fclose(fid_out);    
    end
end





   