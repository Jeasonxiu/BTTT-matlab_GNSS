function [GNGGA] = getGNGGA(filename)
%
% function [GNGGA] = getGNGGA(filename)
%
%   Read the given Logged NMEA file, get GNGGA from NMEA
%   
%   input filename : logged NMEA file
%
%   Example : [GNGGA] = getGNGGA('NMEA.txt')
%
%   coded by Joonseong Gim, Junr 1, 2016
%   GNGGA의 값이 없는 경우에는 무시하고 GNGGA를 추출하여 저장하는 기능 추가
% filename = 'test.txt';
fid=fopen(filename,'r');
% fid = fopen('ubx2_150707_10m_nmea.txt','r');
if fid == -1
    disp('Cannot locate the input file!')
    GNGGA = {};
else
    GNGGA = {};
    fid_out = fopen('GNGGA.txt', 'w');
    GNGGAlist = textscan(fid,'%s');
    stop = 0;
    for i = 1:length(GNGGAlist{1})
        line = cell2mat(GNGGAlist{1}(i));
        if length(line) >= 6
            check = line(1:6);
            if check == '$GNGGA'
                fprintf(fid_out, '%s \n',line);
                GNGGA(i,1) = {line};
            end
        end
    end
    if ~isempty(GNGGA)
        GNGGA = GNGGA(find(~cellfun(@isempty,GNGGA)),1);
        fclose(fid_out);
    end
end


           