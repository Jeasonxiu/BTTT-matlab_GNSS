function [GPGGA] = getGPGGA(filename)
%
% function [GPGGA] = getGPGGA(filename)
%
%   Read the given Logged NMEA file, get GPGGA from NMEA
%
%   input filename : logged NMEA file
%
%   Example : [GPGGA] = getGPGGA('NMEA.txt')
%
%   coded by Joonseong Gim, Jan 21, 2016
%   GPGGA의 값이 없는 경우에는 무시하고 GPGGA를 추출하여 저장하는 기능 추가
% filename = 'test.txt';
fid=fopen(filename,'r');
% fid = fopen('ubx2_150707_10m_nmea.txt','r');
if fid == -1
    disp('Cannot locate the input file!')
    GPGGA = {};
else
    GPGGA = {};
    %     fid_out = fopen('GPGGA.txt', 'w');
    GPGGAlist = textscan(fid,'%s');
    stop = 0;
    cnt=1;
    for i = 1:length(GPGGAlist{1})
%         for i = 27:30
        line = cell2mat(GPGGAlist{1}(i));
        if length(line) >= 6
            check = line(1:6);
            check2 = line(1:5);
            if check == '$GPGGA' %check(1:5) == 'GPGGA'
                %                 fprintf(fid_out, '%s \n',line);
                GPGGA(cnt,1) = {line};
                cnt=cnt+1;
            elseif check == '$GNGGA'
                %                 fprintf(fid_out, '%s \n',line);
                GPGGA(cnt,1) = {line};
                cnt=cnt+1;
            end
            if check2 == 'GPGGA'
                GPGGA(cnt,1) = {strcat('$',line)};
                cnt=cnt+1;
            elseif check2 == 'GNGGA'
                %                 fprintf(fid_out, '%s \n',line);
                
                GPGGA(cnt,1) = {strcat('$',line)};
                cnt=cnt+1;
            end
        end
    end
    if ~isempty(GPGGA)
        GPGGA = GPGGA(find(~cellfun(@isempty,GPGGA)),1);
        %         fclose(fid_out);
    end
end


