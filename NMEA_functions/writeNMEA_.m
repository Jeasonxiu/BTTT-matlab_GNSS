function [NMEAQM,nmeaqm] = writeNMEA_(filename)
    
    % function [NMEAQM] = writeNMEA(filename)
    %
    %   Read the NMEAfile, make a NMEAQM matrix from file
    %
    %   input filename : NMEA file
    %
    %   output NMEAQM : SA by 8 matrix total epoch(gs, lati, longi, alti, prn, el, az, snr)
    %
    %   Example : [NMEAQM] = writeNMEA('jprA014a.txt')
    %
    %   coded by Joonseong Gim, Jan 26, 2016

    
%     clear all;
%     filename = 'DBUU047j_3.ubx';
%     filename = 'jprA063x_opG.txt';
%     filename = 'DBUU049i_16.ubx';
%     filename = 'S7170912.txt';
    %% NMEA유무에 따라 동작
    fid=fopen(filename,'r');
    if fid == -1
        disp('Cannot locate the input file!')
        NMEAQM ={};
    else
        GPRMC = getGPRMC(filename);
        if ~isempty(GPRMC)
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
        else
            yyyy = 0; mm = 0; dd = 0;
            YMD = [yyyy mm dd];
        end
        
        NMEA = getNMEA(filename);
        for i = 1:length(NMEA)
            epoch = NMEA{i,1};
            check = cell2mat(epoch(1));   % GPGGA
            leng = length(epoch);
            if length(check) >= 70
                if leng > 1
                    check2 = cell2mat(epoch(2));   % GPGSA
                    check2 = check2(1:6);
                    if check2 == '$GPGSA'
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
    
    
