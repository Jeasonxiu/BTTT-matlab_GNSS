function [nmea] = writeNMEA(filename,year,month,day)
%
% function [nmea] = writeNMEA(filename,year,month,day))
%
%   Read the NMEAfile, make a NMEAQM matrix from file
%
%   input filename : ubx file(logging with NMEA)
%
%   output NMEAQM : SA by 8 matrix total epoch(gs, lati, longi, alti, prn, el, az, snr)
%
%   Example : [nmea] = writeUBXNMEA('jprA014a.ubx',2017,03,20)
%
%   coded by Joonseong Gim, June 2, 2016
%   Modified by Joonseong Gim, Mar 20, 2017


% clear all;
% filename = 'rover_D_160601_6.ubx';
% filename = '1209_1621_C_GEODE';
% filename = 'GR_170320_jamsil.ubx';
% filename = '20170320.nmea';
% filename = 'ublox_joon.ubx';
% year = 2017; month =03; day = 20;
%% NMEA유무에 따라 동작
fid=fopen(filename,'r');
if fid == -1
    disp('Cannot locate the input file!')
    UBXNMEAQM ={};
else
    [NMEAlist] = NMEALIST(filename);
    [UBXNMEA] = ubxNMEA(NMEAlist);
    cnt=1;
    GSV ={};
    for i=1:length(UBXNMEA)
        epoch = UBXNMEA{i};
        gsa = 0; gsv =0 ;
        for j = 1:length(epoch)
            line = cell2mat(epoch(j));
            if line(2:6) == 'GPGGA' | line(2:6) == 'GNGGA'
                GGA = textscan(line(1:end-4),'%s %s %f %s %f %s %f %f %f %f %s %f %s %f %f','delimiter',',');
                UTC=cell2mat(GGA{1,2});
                if length(UTC) > 5
                    hh = str2double(UTC(1:2));
                    mm = str2double(UTC(3:4));
                    ss = str2double(UTC(5:end));
                    [gw, gs] = date2gwgs(year, month, day, hh, mm, ss);
                end
            elseif line(2:6) == 'GNGSA' | line(2:6) == 'GPGSA'
                gsa = gsa+1;
                GSA(gsa,:) = textscan(line(1:end-3),'%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','delimiter',',');
            elseif line(2:6) == 'GPGSV' | line(2:6) == 'GLGSV'
                gsv = gsv+1;
                GSV(gsv,:) = textscan(line(1:end-3),'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','delimiter',',');
            end
        end
        %% GGA 위치 값이 존재할때에만 행렬을 생성
        if ~isnan(cell2mat(GGA(3)))
            %% Latitude 생성
            degLat = fix(cell2mat(GGA(3))/100);
            minLat = cell2mat(GGA(3))-degLat*100;
            la = degLat + minLat/60;        % Latitude
            %% Longitude 생성
            degLon = fix(cell2mat(GGA(5))/100);
            minLon = cell2mat(GGA(5))-degLon*100;
            lon = degLon + minLon/60;        % Latitude
            %% Hegiht 생성
            ht = cell2mat(GGA(10)) + cell2mat(GGA(12));
            nmea(cnt,1:4) = [round(gs), la, lon, ht];
            if ~isempty(GSA)
                for k = 1:gsa
                    for kk =4:15
                        prn = cell2mat(GSA(k,kk));
                        if ~isnan(prn)
                            nmea(cnt,1:5) = [round(gs), la, lon, ht, prn];
                            cnt=cnt+1;
                        end
                    end
                end
            else
                cnt=cnt+1;
            end
            if ~isempty(GSV)
                for k=1:gsv
                    for kk=5:4:17
                        prn = cell2mat(GSV(k,kk));
                        az = [cell2mat(GSV(k,kk+1))]; el = [cell2mat(GSV(k,kk+2))]; snr = [cell2mat(GSV(k,kk+3))];
                        if isempty(az) | isnan(az)
                            az = 0; end
                        if isempty(el) | isnan(el)
                            el = 0; end
                        if isempty(snr) | isnan(snr)
                            snr = 0; end
                        if ~isnan(prn)
                            same = max(find(nmea(:,5) == prn));
                            if ~isempty(same)
                                nmea(same,6:8) = [az, el, snr];
                            end
                        end
                    end
                end
            end
        end
        
    end
      
end
