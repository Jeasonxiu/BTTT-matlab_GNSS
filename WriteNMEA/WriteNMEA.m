function [nmea] = WriteNMEA(filename,year,month,day)
%
% function [nmea] = WriteNMEA(filename,year,month,day))
%
%   Read the NMEAfile, make a NMEAQM matrix from file
%
%   input filename : ubx or any text file(include NMEA message)
%
%   output NMEAQM : N by 10 matrix
%                  (gs, latitude, longitude, height, Status, Number of Sats, prn, el, az, snr)
%                 prn : GPS 100, BDS 200, GLO 300, QZSS 500, other 600
%
%   Example : [nmea] = writeNMEA('jprA014a.ubx',2017,03,20)
%
%   coded by Joonseong Gim, June 2, 2016
%   Modified by Joonseong Gim, Mar 20, 2017
%   Beidou Az, El, SNR 추가
%   Modified by Joonseong Gim, Mar 29, 2017
%   QZSS Az, El, SNR 추가
%   Modified by Joonseong Gim, Jan 02, 2018


%% test 용 파일 및 파라미터들
% clear all;
% filename = '362NMEA.txt';
% year = 2017; month = 12; day = 28;

%% NMEA유무에 따라 동작
fid=fopen(filename,'r');
if fid == -1
    disp('No file!')
    UBXNMEAQM ={};
else
    [NMEAlist] = NMEALIST(filename);        % 에폭별 NMEA 리스트업
    [UBXNMEA] = ubxNMEA(NMEAlist);          % 에폭별 NMEA 리스트 셀어레이
    cnt=1;
    GSA ={};GSV ={}; prn = 0;
    for i=1:length(UBXNMEA)
        epoch = UBXNMEA{i};
        la =0; gsa = 0; gsv =0 ; GGA{1,3} = nan; RMC{1,4} = nan; HT = {};
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
                FQ = cell2mat(GGA(7));
                NSV = cell2mat(GGA(8));
            elseif line(2:6) == 'GPRMC' | line(2:6) == 'GNRMC'
                RMC = textscan(line(1:end-4),'%s %s %s %f %s %f %s %f %f %s %s %s','delimiter',',');
                UTC = cell2mat(RMC{1,2});
                if length(UTC) > 5
                    hh = str2double(UTC(1:2));
                    mm = str2double(UTC(3:4));
                    ss = str2double(UTC(5:end));
                    [gw, gs] = date2gwgs(year, month, day, hh, mm, ss);
                end
            elseif line(2:6) == 'GNGSA' | line(2:6) == 'GPGSA' | line(2:6) == 'BDGSA' |...
                    line(2:6) == 'QZGSA'
                gsa = gsa+1;
                GSA(gsa,:) = textscan(line(1:end-3),...
                    '%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','delimiter',',');
                %% GPGSA와 GNGPA의 GPS 사용 정보가 같을때 GNGSA만 사용하도록 설정
                if gsa > 1
                    if cell2mat(GSA(gsa,3:6)) == cell2mat(GSA(gsa-1,3:6))
                        gsa = gsa - 1; end; end
            elseif line(2:6) == 'GPGSV' | line(2:6) == 'GLGSV' | line(2:6) == 'BDGSV' |...
                    line(2:6) == 'GBGSV' | line(2:6) == 'QZGSV'
                gsv = gsv+1;
                
                if length(line) < 80
                    GSV(gsv,:) = textscan(line(1:end-3),...
                        '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','delimiter',',');
                elseif length(line)> 80
                    GSV(gsv,:) = textscan(line(1:70-3),...
                        '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','delimiter',',');
                end
                
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
            HT = [ht];
        elseif ~isnan(cell2mat(RMC(4)))
            %% Latitude 생성
            degLat = fix(cell2mat(RMC(4))/100);
            minLat = cell2mat(RMC(4))-degLat*100;
            la = degLat + minLat/60;        % Latitude
            %% Longitude 생성
            degLon = fix(cell2mat(RMC(6))/100);
            minLon = cell2mat(RMC(6))-degLon*100;
            lon = degLon + minLon/60;        % Latitude
            %% Hegiht 생성
            if isnan(cell2mat(GGA(3))) && length(GGA) < 4
                ht = 0; FQ = 0; NSV =0;
            end; end
        if la ~= 0
            nmea(cnt,1:6) = [round(gs), la, lon, ht, FQ, NSV];
            
            if ~isempty(GSA)
                for k = 1:gsa
                    for kk =4:15
                        prn = cell2mat(GSA(k,kk));
                        System = cell2mat(GSA{k,1});
                        if System(2:6) == 'QZGSA'
                            prn = prn + 500;
                        end
                        
                        if ~isnan(prn)
                            if prn > 140 & prn < 200
                                prn = prn + 60;
                            elseif prn > 64 & prn < 96
                                prn = prn + 236;
                            elseif prn > 37 & prn < 60
                                prn = prn + 263;
                            elseif prn < 100
                                prn = prn + 100;
                            elseif prn > 200
                                prn = prn;
                            else
                                prn = prn + 600;
                            end
                            if prn < 600
                                nmea(cnt,1:7) = [round(gs), la, lon, ht, FQ, NSV, prn];
                                cnt=cnt+1; end
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
                        System = cell2mat(GSV{k,1});
                        if System(2:6) == 'QZGSV'
                            prn = prn + 500;
                        end
                        if prn > 140 & prn < 200
                            prn = prn + 60;
                        elseif prn > 64 & prn < 96
                            prn = prn + 236;
                        elseif prn > 37 & prn < 60
                            prn = prn + 263;
                        elseif prn > 100
                            prn = prn;
                        elseif prn > 200
                            prn = prn;
                        else
                            prn = prn + 100;
                        end
                        az = [cell2mat(GSV(k,kk+1))]; el = [cell2mat(GSV(k,kk+2))]; snr = [cell2mat(GSV(k,kk+3))];
                        if isempty(az) | isnan(az)
                            az = 0; end
                        if isempty(el) | isnan(el)
                            el = 0; end
                        if isempty(snr) | isnan(snr)
                            snr = 0; end
                        if ~isnan(prn) & length(nmea(1,:)) > 4
                            same = max(find(nmea(:,7) == prn));
                            if ~isempty(same)
                                nmea(same,8:10) = [az, el, snr];
                            end
                        end
                    end
                end
            end
        end
        
    end
end
if length(find(nmea(:,5)==0)) ~= 0
    nmea = nmea(max(find(nmea(:,5)==0))+1:end,:);
else
end