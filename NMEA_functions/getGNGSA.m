function [GNGSA] = getGNGSA(filename)
%
% function [GNGSA] = getGNGSA(filename)
%
%   Read the given Logged NMEA file, get GPGSV from NMEA
%   
%   input filename : logged NMEA file
%
%   Example : [] = getGNGSA('NMEA.txt')
%
%   coded by Joonseong Gim, June 1, 2016
%   
% GPS DOP and active satellites
% 
% $GNGSA,A,3,01,07,08,27,11,16,30,,,,,,2.23,0.90,2.04*1F
% $GNGSA,A,3,84,80,83,74,,,,,,,,,2.23,0.90,2.04*1C
% 
% 



fid = fopen(filename,'r');
% fid = fopen('ubx2_150707_10m_nmea.txt');

if fid == -1
    disp('Cannot locate the input file!')
    GNGSA = {};
else
    GNGSA = {};
    fid_out = fopen('GPGSV.txt', 'w');
    GNGSAlist = textscan(fid,'%s');
    
    for j = 1:length(GNGSAlist{1})
        line = cell2mat(GNGSAlist{1}(j));
        if length(line) >= 6
            nmea = line(1:6);
            if nmea == '$GNGSA'
                fprintf(fid_out, '%s \n',line);
                coloum(j,1) = j;
                GNGSA(j,1) = {line};
            end
        end
    end
    if ~isempty(GNGSA)
        GNGSA = GNGSA(find(~cellfun(@isempty,GNGSA)),1);
        fclose(fid_out);    
    end
end

% for i = 1:length(GNGSA)
%     line2 = cell2mat(GNGSA{1}(12:13));
%     if line2 == ',,'
%     else
%         prn = str2num(line2);
%     end
%     if ~isempty(prn)
%         if prn < 33
%             GSA{i,1}
%         
%     
% 
% 
% 
%    