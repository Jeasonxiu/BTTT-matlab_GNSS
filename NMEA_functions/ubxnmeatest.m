% function [NMEAQM,nmeaqm] = writeUBXNMEA(filename)
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


clear all;
% filename = 'rover_D_160601_6.ubx';
filename = 'DRUA1531.ubx';

%% NMEA유무에 따라 동작
fid=fopen(filename,'r');
if fid == -1
    disp('Cannot locate the input file!')
    UBXNMEAQM ={};
else
    [NMEAlist] = NMEALIST(filename);
    [UBXNMEA] = ubxNMEA(filename);
    %% gs를 만들기위해 RMC 유무 확인
    GPRMC = ubxGPRMC(NMEAlist);
    GNRMC = ubxGNRMC(NMEAlist);
    GPGGA = ubxGPGGA(NMEAlist);
    GLGSV = ubxGLGSV(NMEAlist);
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
    %% 단독항법 NMEA인지 다중항법인지 판단하기 위해 GPGGA 유무 확인
    if ~isempty(GLGSV)
        GGA = 'GPGGA';
        GSA = 'GPGSA';
        sys = 1;
    else
        GGA = 'GNGGA';
        GSA = 'GNGSA';
        sys = 2;
    end
    
    switch sys
        case 1
            for i = 1:length(UBXNMEA(:,1))
                epoch = UBXNMEA{i,1};
                check = cell2mat(epoch(1));   % GNGGA
                leng = length(epoch);
                if length(check) >= 70
                    if leng > 1
                        check2 = cell2mat(epoch(2));   % GNGSA
                        check2 = check2(2:6);
                         if check2 == GSA
                        switch leng
                            case 6
                                gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
                                gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                gpgsv = epoch(3:6,1);
                                nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
                                nmeaqm(i,1) = {nmea};
                            case 5
                                gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
                                gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                gpgsv = epoch(3:5,1);
                                nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
                                nmeaqm(i,1) = {nmea};
                            case 4
                                gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
                                gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                gpgsv = epoch(3:4,1);
                                nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
                                nmeaqm(i,1) = {nmea};
                            case 3
                                gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
                                gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                gpgsv = epoch(3,1);
                                nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
                                nmeaqm(i,1) = {nmea};
                            case 2
                                gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
                                gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                gpgsv = [];
                                nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
                                nmeaqm(i,1) = {nmea};
                        end
                    elseif check2 == '$GPGSV'
                        switch leng
                            case 5
                                gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
                                gpgsa = [];
                                gpgsv = epoch(2:5,1);
                                nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
                                nmeaqm(i,1) = {nmea};
                            case 4
                                gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
                                gpgsa = [];
                                gpgsv = epoch(2:4,1);
                                nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
                                nmeaqm(i,1) = {nmea};
                            case 3
                                gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
                                gpgsa = [];
                                gpgsv = epoch(2:3,1);
                                nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
                                nmeaqm(i,1) = {nmea};
                            case 2
                                gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
                                gpgsa = [];
                                gpgsv = epoch(2,1);
                                nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
                                nmeaqm(i,1) = {nmea};
                        end
                    end
                else
                    gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
                    gpgsa = [];
                    gpgsv = [];
                    nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
                    nmeaqm(i,1) = {nmea};
                end
            end
        end
        case 2
            for i = 1:length(UBXNMEA(:,1))
                epoch = UBXNMEA{i,1};
                check = cell2mat(epoch(1));   % GNGGA
                leng = length(epoch);
                if length(check) >= 70
                    if leng > 1
                        check2 = cell2mat(epoch(2));   % GNGSA
                        check2 = check2(2:6);
                        if check2 == GSA
                            gp = cell2mat(epoch(4));   % GPGSV
                            numgp = str2num(gp(8));
                            gl = cell2mat(epoch(length(epoch)));    % GLGSV
                            numgl = str2num(gl(8)); 
                            num_line = 3 + numgp + numgl;
                            switch num_line
                                case 11
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                    glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
                                    gpgsv = epoch(4:3+numgp,1);
                                    glgsv = epoch(3+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 10
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                    glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
                                    gpgsv = epoch(4:3+numgp,1);
                                    glgsv = epoch(3+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 9
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                    glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
                                    gpgsv = epoch(4:3+numgp,1);
                                    glgsv = epoch(3+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 8
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                    glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
                                    gpgsv = epoch(4:3+numgp,1);
                                    glgsv = epoch(3+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 7
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                    glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
                                    gpgsv = epoch(4:3+numgp,1);
                                    glgsv = epoch(3+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 6
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                    glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
                                    gpgsv = epoch(4:3+numgp,1);
                                    glgsv = epoch(3+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 5
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                    glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
                                    gpgsv = epoch(4:3+numgp,1);
                                    glgsv = epoch(3+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 4
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                    glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
                                    gpgsv = epoch(4:3+numgp,1);
                                    glgsv = epoch(3+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 3
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                    glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
                                    gpgsv = epoch(4:3+numgp,1);
                                    glgsv = epoch(3+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 2
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                    glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
                                    gpgsv = epoch(4:3+numgp,1);
                                    glgsv = epoch(3+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 1
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
                                    glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
                                    gpgsv = epoch(4:3+numgp,1);
                                    glgsv = epoch(3+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                            end
                        elseif check2 == 'GPGSV'
                            switch num_line
                                case 11
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = [];
                                    glgsa = [];
                                    gpgsv = epoch(2:1+numgp,1);
                                    glgsv = epoch(1+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 10
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = [];
                                    glgsa = [];
                                    gpgsv = epoch(2:1+numgp,1);
                                    glgsv = epoch(1+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 9
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = [];
                                    glgsa = [];
                                    gpgsv = epoch(2:1+numgp,1);
                                    glgsv = epoch(1+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 8
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = [];
                                    glgsa = [];
                                    gpgsv = epoch(2:1+numgp,1);
                                    glgsv = epoch(1+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 7
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = [];
                                    glgsa = [];
                                    gpgsv = epoch(2:1+numgp,1);
                                    glgsv = epoch(1+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 6
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = [];
                                    glgsa = [];
                                    gpgsv = epoch(2:1+numgp,1);
                                    glgsv = epoch(1+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 5
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = [];
                                    glgsa = [];
                                    gpgsv = epoch(2:1+numgp,1);
                                    glgsv = epoch(1+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                                case 4
                                    gga = strsplit(cell2mat(epoch(1)),{',','*'});
                                    gpgsa = [];
                                    glgsa = [];
                                    gpgsv = epoch(2:1+numgp,1);
                                    glgsv = epoch(1+numgp+1:leng,1);
                                    nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                                    nmeaqm(i,1) = {nmea};
                            end
                        end
                    else
                        gga = strsplit(cell2mat(epoch(1)),{',','*'});
                        gpgsa = [];
                        glgsa = [];
                        gpgsv = [];
                        glgsv = [];
                        nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
                        nmeaqm(i,1) = {nmea};
                    end
                end
            end
    end
    
    nmeaqm = nmeaqm(find(~cellfun(@isempty,nmeaqm)),1);
    
    i=1;
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
end
