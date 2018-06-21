function [NMEAQM,nmeaqm] = writeNMEA2(filename)
%
% function [NMEAQM] = writeNMEA(filename)
%
%   Read the NMEAfile, make a NMEAQM matrix from file
%
%   input filename : ubx file(logging with NMEA)
%
%   output NMEAQM : SA by 8 matrix total epoch(gs, lati, longi, alti, prn, el, az, snr)
%
%   Example : [NMEAQM] = writeUBXNMEA('jprA014a.ubx')
%
%   coded by Joonseong Gim, June 2, 2016


% clear all;
% filename = 'rover_D_160601_6.ubx';
% filename = 'jprA214f_s7.txt';
% filename = 'jprA014a_tab.txt';
% filename = 'GR_170320_jamsil.ubx';
%% NMEA유무에 따라 동작
fid=fopen(filename,'r');
if fid == -1
    disp('Cannot locate the input file!')
    UBXNMEAQM ={};
else
    [NMEAlist] = NMEALIST(filename);
    [UBXNMEA] = ubxNMEA(NMEAlist);
end
GPGGA = {}; GNGGA = {}; GPGSA = {}; QZGSA = {}; GPGSV = {}; GLGSV = {}; GNGSA = {}; GPRMC ={}; GNRMC ={};
for i = 1:length(UBXNMEA)
    EPOCH = UBXNMEA{i,1};
    cnt = 1; cnt2 = 1; cnt3 = 1;
    for j = 1:length(EPOCH)
        line = cell2mat(EPOCH(j));
        if line(2:6) == 'GPGGA'
            split = strsplit(line,',');
            if length(split) > 3
                used_sat = str2num(cell2mat(split(8)));
                GPGGA{i,1} = line;
            end
        end
        if line(2:6) == 'GNGGA'
            split = strsplit(line,',');
            if length(split) > 3
                used_sat = str2num(cell2mat(split(8)));
                GNGGA{i,1} = line;
            end
        end
        if line(2:6) == 'GPGSA'
            split = strsplit(line,',');
            if length(split) > 6
                gps_num_sat = length(split(4:length(split)-3));
                GPGSA{i,1} = line;
            end
        end
        if line(2:6) == 'QZGSA'
            split = strsplit(line,',');
            if length(split) > 6
                qz_num_sat = length(split(4:length(split)-3));
                QZGSA{i,1} = line;
            end
        end
        if line(2:6) == 'GPGSV'
            gpgsv{cnt,1} = line;
            GPGSV{i,1} = gpgsv;
            cnt = cnt + 1;
        end
        if line(2:6) == 'GLGSV'
            glgsv{cnt2,1} = line;
            GLGSV{i,1} = glgsv;
            cnt2 = cnt2 + 1;
        end
        if line(2:6) == 'GNGSA'
            gngsa{cnt3,1} = line;
            GNGSA{i,1} = gngsa;
            cnt3 = cnt3 + 1;
        end
        if line(2:6) == 'GPRMC'
            GPRMC{i,1} = line;
        end
        if line(2:6) == 'GNRMC'
            GNRMC{i,1} = line;
        end
    end
    
    if used_sat == gps_num_sat
        gga = strsplit(GPGGA{i},{',','*'});
        if ~isempty(GPRMC{i})
            [YMD] = RMC2YMD(GPRMC{i});
        else
            YMD = [0 0 0];
        end
        if ~isempty(GPGSA{i}) && ~isempty(GPGSV{i})
            gsa = strsplit(GPGSA{i},{',','*'});
            gsv = GPGSV{i};
            nmea = mknmea(YMD, gga, gsa, gsv);
            nmeaqm(i,1) = {nmea};
        elseif ~isempty(GPGSA{i})
            gsa = strsplit(GPGSA{i},{',','*'});
            gsv = [];
            nmea = mknmea(YMD, gga, gsa, gsv);
            nmeaqm(i,1) = {nmea};
        elseif ~isempty(GPGSV{i})
            gsa = [];
            gsv = GPGSV{i};
            nmea = mknmea(YMD, gga, gsa, gsv);
            nmeaqm(i,1) = {nmea};
        else
            gsa = [];
            gsv = [];
            nmea = mknmea(YMD, gga, gsa, gsv);
            nmeaqm(i,1) = {nmea};
        end
    else
        if ~isempty(GNGGA)
            gga = strsplit(GNGGA{i},{',','*'});
        else
            gga = strsplit(GPGGA{i},{',','*'});
        end
        if ~isempty(GNRMC)
            [YMD] = RMC2YMD(GNRMC{i});
        elseif ~isempty(GPRMC{i})
            [YMD] = RMC2YMD(GPRMC{i});
        else
            YMD = [0 0 0];
        end
        
        if ~isempty(GNGSA{i}) && ~isempty(GPGSV{i}) && ~isempty(GLGSV{i})
            gngsa = GNGSA{i};
            gpgsv = GPGSV{i};
            glgsv = GLGSV{i};
            nmea = mknmea2(YMD, gga, gngsa, gpgsv, glgsv);
            nmeaqm(i,1) = {nmea};
        elseif ~isempty(GPGSA{i})
            gsa = strsplit(GPGSA{i},{',','*'});
            gsv = [];
            nmea = mknmea2(YMD, gga, gngsa, gpgsv, glgsv);
            nmeaqm(i,1) = {nmea};
        elseif ~isempty(GPGSV{i})
            gsa = [];
            gsv = GPGSV{i};
            nmea = mknmea2(YMD, gga, gngsa, gpgsv, glgsv);
            nmeaqm(i,1) = {nmea};
        else
            gsa = [];
            gsv = [];
            nmea = mknmea2(YMD, gga, gngsa, gpgsv, glgsv);
            nmeaqm(i,1) = {nmea};
        end
        gps_view = strsplit(cell2mat(GPGSV{i}(1)),',');
        gpsview(i,1) = str2num(cell2mat(gps_view(4)));
        glo_view = strsplit(cell2mat(GLGSV{i}(1)),',');
        gloview(i,1) = str2num(cell2mat(glo_view(4)));
        qz_view = strsplit(QZGSA{i},',');
        qzview(i,1) = str2num(cell2mat(qz_view(4)));
    end
    
end
for k =1:length(GNGSA)
    GSA = GNGSA{k};
    gngsa1 = strsplit(cell2mat(GSA(1)),{',','*'}); gngsa2 = strsplit(cell2mat(GSA(2)),{',','*'});
    GSA_num(k,:) = [length(gngsa1(4:length(gngsa1)-4)), length(gngsa2(4:length(gngsa2)-4))];
end
i = 1; j = 1;
for j = 1:length(nmeaqm)
    if length(nmeaqm{j}) >3
        for jj = 1:length(nmeaqm{j}(:,1))
            NMEAQM(i,1:length(nmeaqm{j}(1,:))) = nmeaqm{j}(jj,1:length(nmeaqm{j}(1,:)));
            i = i + 1;
        end
    else
        NMEAQM(j,1:3) = nmeaqm{j};
    end
end