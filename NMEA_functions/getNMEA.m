function [NMEA] = getNMEA(filename)
%
% function [GPGSA] = getGPGSA(filename)
%
%   Read the given Logged NMEA file, get GPGGA, GPGSA, GPGSV from logged
%   file
%
%   input filename : logged NMEA file
%   output NMEA : GPGGA, GPGSA, GPGSV cell array / epoch
%
%   Example : [] = getNMEA('NMEA.txt')
%
%   coded by Joonseong Gim, Jan 25, 2016
%


% clear all;
% filename = 'DBUU047j_3.ubx';
% filename = 'jprA014a.txt';
% filename = '(b)20160112120211.txt';
%% NMEA파일로 부터 GPGGA, GPGSV, GPGSA, GPRMC를 추출하는 과정
fid=fopen(filename,'r');
if fid == -1
    disp('Cannot locate the input file!')
    NMEA = {};
else
    GPGSV = getGPGSV(filename);
    GPGGA = getGPGGA(filename);
    GPGSA = getGPGSA(filename);
    if ~isempty(GPGSA)
        if length(GPGSA) >= length(GPGGA)
            a = {}; b = {}; % NMEA matrix 의 GGA, GSA 를 위한 빈 Cell 생성
            if ~isempty(GPGSV)
                %% GPGSV에서 Seq 넘버와 Msg 넘버를 가지고 한 Epoch당 GPGSV Cell array 생성
                for i = 1: length(GPGSV)
                    line = cell2mat(GPGSV(i,1));
                    index = findstr(line,',');
                    ToNo = str2num(line(index(1)+1:index(2)-1));
                    msgNo = str2num(line(index(2)+1:index(3)-1));
                    if ToNo == 4 && msgNo == 1
                        line2 = cell2mat(GPGSV(i+1,1));
                        line3 = cell2mat(GPGSV(i+2,1));
                        line4 = cell2mat(GPGSV(i+3,1));
                        epoch{i,1} = {a; b; line; line2; line3; line4};
                    elseif ToNo == 3 && msgNo == 1
                        line2 = cell2mat(GPGSV(i+1,1));
                        line3 = cell2mat(GPGSV(i+2,1));
                        epoch{i,1} = {a; b; line; line2; line3};
                    elseif ToNo == 2 && msgNo == 1
                        line2 = cell2mat(GPGSV(i+1,1));
                        epoch{i,1} = {a; b; line; line2};
                    else
                    end
                end
                GPGSV = epoch(find(~cellfun(@isempty,epoch)),1);
                NMEA = GPGSV;
                for j = 1: length(GPGGA)
                    NMEA{j,1}(1) = {cell2mat(GPGGA(j,1))};
                    NMEA{j,1}(2) = {cell2mat(GPGSA(j,1))};
                end
            else
                for j = 1: length(GPGGA)
                    NMEA{j,1}(1,1) = {cell2mat(GPGGA(j,1))};
                    NMEA{j,1}(2,1) = {cell2mat(GPGSA(j,1))};
                end
            end
            
            for j = 1: length(GPGGA)
                NMEA{j,1}(1,1) = {cell2mat(GPGGA(j,1))};
                NMEA{j,1}(2,1) = {cell2mat(GPGSA(j,1))};
            end
            
        elseif ~isempty(GPGSV)
            a = {}; % NMEA matrix 의 GGA 를 위한 빈 Cell 생성
            for i = 1: length(GPGSV)
                line = cell2mat(GPGSV(i,1));
                index = findstr(line,',');
                ToNo = str2num(line(index(1)+1:index(2)-1));
                msgNo = str2num(line(index(2)+1:index(3)-1));
                if ToNo == 4 && msgNo == 1
                    line2 = cell2mat(GPGSV(i+1,1));
                    line3 = cell2mat(GPGSV(i+2,1));
                    line4 = cell2mat(GPGSV(i+3,1));
                    epoch{i,1} = {a; line; line2; line3; line4};
                elseif ToNo == 3 && msgNo == 1
                    line2 = cell2mat(GPGSV(i+1,1));
                    line3 = cell2mat(GPGSV(i+2,1));
                    epoch{i,1} = {a; line; line2; line3};
                elseif ToNo == 2 && msgNo == 1
                    line2 = cell2mat(GPGSV(i+1,1));
                    epoch{i,1} = {a; line; line2; line3};
                else
                end
            end
            GPGSV = epoch(find(~cellfun(@isempty,epoch)),1);
            NMEA = GPGSV;
            for j = 1: length(GPGGA)
                NMEA{j,1}(1) = {cell2mat(GPGGA(j,1))};
            end
        else
            for j = 1: length(GPGGA)
                NMEA{j,1}(1) = {cell2mat(GPGGA(j,1))};
            end
        end
    end
end



