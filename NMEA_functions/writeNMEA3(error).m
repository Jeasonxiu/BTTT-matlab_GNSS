function [NMEAQM,nmeaqm] = writeNMEA3(filename)
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
% filename = 'jprB050h_tab.txt';
% filename = 'jprA032f_s6.txt';
% filename = 'GR_170320_jamsil.ubx';
% filename = 'PPS1_170407.NMEA';
% filename = 'test4.cnb';
%% NMEA유무에 따라 동작
fid=fopen(filename,'r');
if fid == -1
    disp('Cannot locate the input file!')
    UBXNMEAQM ={};
else
    [NMEAlist] = NMEALIST(filename);
    [UBXNMEA] = ubxNMEA(NMEAlist);
    %% gs를 만들기위해 RMC 유무 확인
    GPRMC = ubxGPRMC(NMEAlist);
    GNRMC = ubxGNRMC(NMEAlist);
    GPGGA = ubxGPGGA(NMEAlist);
    GNGSA = ubxGNGSA(NMEAlist);
    if ~isempty(GPRMC)
        RMC = 1;
        GPRMC = GPRMC;
    elseif ~isempty(GNRMC)
        RMC = 2;
        GPRMC = GNRMC;
    else
        RMC = 3;
    end
    %% RMC가 존재할 경우 yyyy, mm, dd 생성
    if RMC < 3
        for i = 1: length(GPRMC)
            line = cell2mat(GPRMC(i,1));
            if length(line) >= 60
                break
            end
        end
        if length(line) >= 60
            index = findstr(line,',');
            yymmdd = line(index(9)+1:index(10)-1);
            if str2num(yymmdd(5:6)) >= 80
                yyyy = 1900 + str2num(yymmdd(5:6));
            else
                yyyy = 2000 + str2num(yymmdd(5:6));
            end
            mm = str2num(yymmdd(3:4));
            dd = str2num(yymmdd(1:2));
            YMD = [yyyy mm dd];
        else
            yyyy = 0; mm = 0; dd = 0;
            YMD = [yyyy mm dd];
        end
    end
    for i = 1:length(UBXNMEA(:,1))
        epoch = UBXNMEA{i,1};
        check = cell2mat(epoch(1));   % GNGGA
        nmea_list = {};
        
        for j = 1:length(epoch(:,1))
            nmea_list{j,1} = epoch{j}(1,2:6);
        end
        nmea_list = cell2mat(nmea_list);
        numofgpgsv = strmatch('GPGSV',nmea_list);
        numofglgsv = strmatch('GLGSV',nmea_list);
        numofgpgsa = strmatch('GPGSA',nmea_list);
        numofgngsa = strmatch('GNGSA',nmea_list);
        if length(numofgpgsv) == 4
            gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
            gpgsa = strsplit(cell2mat(epoch(numofgpgsa)),{',','*'});
            gpgsv = epoch(min(numofgpgsv):max(numofgpgsv),1);
            nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
            gpgsvqm(i,1) = {nmea};
        elseif length(numofgpgsv) == 3
            gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
            gpgsa = strsplit(cell2mat(epoch(numofgpgsa)),{',','*'});
            gpgsv = epoch(min(numofgpgsv):max(numofgpgsv),1);
            nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
            gpgsvqm(i,1) = {nmea};
        elseif length(numofgpgsv) == 2
            gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
            gpgsa = strsplit(cell2mat(epoch(numofgpgsa)),{',','*'});
            gpgsv = epoch(min(numofgpgsv):max(numofgpgsv),1);
            nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
            gpgsvqm(i,1) = {nmea};
        else
            gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
            gpgsa = strsplit(cell2mat(epoch(numofgpgsa)),{',','*'});
            gpgsv = epoch(numofgpgsv,1);
            nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
            gpgsvqm(i,1) = {nmea};
        end
        
        if length(numofglgsv) == 4
            gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
            gngsa = strsplit(cell2mat(epoch(max(numofgngsa))),{',','*'});
            glgsv = epoch(min(numofglgsv):max(numofglgsv),1);
            nmea = mknmea(YMD, gpgga, gngsa, glgsv);
            glgsvqm(i,1) = {nmea};
        elseif length(numofglgsv) == 3
            gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
            gngsa = strsplit(cell2mat(epoch(max(numofgngsa))),{',','*'});
            glgsv = epoch(min(numofglgsv):max(numofglgsv),1);
            nmea = mknmea(YMD, gpgga, gngsa, glgsv);
            glgsvqm(i,1) = {nmea};
        elseif length(numofglgsv) == 2
            gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
            gngsa = strsplit(cell2mat(epoch(max(numofgngsa))),{',','*'});
            glgsv = epoch(min(numofglgsv):max(numofglgsv),1);
            nmea = mknmea(YMD, gpgga, gngsa, glgsv);
            glgsvqm(i,1) = {nmea};
        else
            gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
            gngsa = strsplit(cell2mat(epoch(max(numofgngsa))),{',','*'});
            glgsv = epoch(numofglgsv,1);
            nmea = mknmea(YMD, gpgga, gngsa, glgsv);
            glgsvqm(i,1) = {nmea};
        end
        nmeaqm{i,1} = [cell2mat(gpgsvqm(i,1)); cell2mat(glgsvqm(i,1))];
        
        i=1;
        for k = 1:length(nmeaqm)
            if length(nmeaqm{k}) >3
                for jj = 1:length(nmeaqm{k}(:,1))
                    NMEAQM(i,1:length(nmeaqm{k}(1,:))) = nmeaqm{k}(jj,1:length(nmeaqm{k}(1,:)));
                    i = i + 1;
                end
            else
                NMEAQM(k,1:3) = nmeaqm{k};
            end
        end
        %         if length(check) >= 70
        %             check2 = cell2mat(epoch(2));
        %             check2 = check2(2:6);
        %             if check2 == 'GPGSV'
        
    end
    
end
%
