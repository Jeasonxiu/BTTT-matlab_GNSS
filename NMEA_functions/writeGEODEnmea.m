% function [NMEAQM] = writeGEODEnmea(filename,yy,mo,dd)

% function [NMEAQM] = writeGEODEnmea(filename,yy,mo,dd)

%   Read the NMEAfile, make a NMEAQM matrix from file
%
%   input filename : GEODE file(logging with NMEA)
%
%   output NMEAQM : n X 13 matrix total epoch(gs, x, y, z, lati, longi, height, number of sats, nmea Quality,prn, el, az, snr)
%
%   Example : [NMEAQM] = writeGEODEnmea('GEODE file',2016,12,09)

%   coded by Joonseong Gim, DEC 12, 2016


clear all;
filename = '1209_1621_C_GEODE';
% % filename = 'jprB014a_tab.txt';
yy= 2016;mo=12;dd=09;

% NMEA유무에 따라 동작
fid=fopen(filename,'r');
if fid == -1
    disp('Cannot locate the input file!')
    UBXNMEAQM ={};
else
    %     [NMEAlist] = NMEALIST(filename);
    [UBXNMEA] = ubxNMEA(filename);
    % gs를 만들기위해 RMC 유무 확인
    k=1;
    for i = 1:length(UBXNMEA(:,1))
        epoch = UBXNMEA{i,1};
        GPRMC = NMEA_extract(epoch,'$GPRMC');
        GNRMC = NMEA_extract(epoch,'$GNRMC');
        GPGGA = NMEA_extract(epoch,'$GPGGA');
        GNGGA = NMEA_extract(epoch,'$GNGGA');
        GPGSA = NMEA_extract(epoch,'$GPGSA');
        GNGSA = NMEA_extract(epoch,'$GNGSA');
        GPGSV = NMEA_extract(epoch,'$GPGSV');
        GLGSV = NMEA_extract(epoch,'$GLGSV');
        
        %% GNGGA 존재 여부 확인
        if  ~isempty(GNGGA)
            GGA = 'GNGGA';
            GPGGA = cell2mat(GNGGA);
            %             gpgga = strsplit(cell2mat(GNGGA),{',','*'});
        else
            GGA = 'GPGGA';
            GPGGA = cell2mat(GPGGA);
            %             gpgga = strsplit(cell2mat(GPGGA),{',','*'});
        end
        %% 단독항법 NMEA인지 다중항법인지 판단하기 위해 GPGGA 유무 확인
        if ~isempty(GLGSV)
            GSA = 'GNGSA';
            sys = 1;
        else
            GSA = 'GPGSA';
            sys = 2;
        end
        switch sys
            case 1
                
                if length(GPGGA) > 70
                    [hh,mm,ss,x,y,z,la,lo,qi,nSats,ht] = NEWreadGGA3(GPGGA) ;
                    [gw, gs] = date2gwgs(yy, mo, dd, hh, mm, ss);
                    [GPGSVmat] = GSVmat2(GPGSV);
                    [GLGSVmat] = GSVmat2(GLGSV);
                    [GNGSAmat] = GSAmat2(GNGSA);
                    GSV = [GPGSVmat;GLGSVmat];
                    for j=1:length(GNGSAmat(find(GNGSAmat(:,1) < 100)))
                        nmeaqm(j,1:9) = [round(gs), x,y,z, la,lo, ht, nSats, qi];
                        prn = GNGSAmat(j);
                        nmeaqm(j,10:13) = GSV(find(GSV(:,1) == prn),:);
                        NMEAQM(k,:) = nmeaqm(j,:);
                        k = k + 1;
                    end
                end
            case 2
                if length(GPGGA) > 70
                    [hh,mm,ss,x,y,z,la,lo,qi,nSats,ht] = NEWreadGGA3(GPGGA) ;
                    [gw, gs] = date2gwgs(yy, mo, dd, hh, mm, ss);
                    [GPGSVmat] = GSVmat2(GPGSV);
                    [GPGSAmat] = GSAmat2(GPGSA);
                    GSV = [GPGSVmat];
                    for j=1:length(GPGSAmat(find(GPGSAmat(:,1) < 100)))
                        nmeaqm(j,1:9) = [round(gs), x,y,z, la,lo, ht, nSats, qi];
                        prn = GNGSAmat(j);
                        nmeaqm(j,10:13) = GSV(find(GSV(:,1) == prn),:);
                        NMEAQM(k,:) = nmeaqm(j,:);
                        k = k + 1;
                    end
                end
        end
    end
end

%     switch sys
%         case 1
%             for i = 1:length(UBXNMEA(:,1))
%                 epoch = UBXNMEA{i,1};
%                 check = cell2mat(epoch(1));   % GNGGA
%                 leng = length(epoch);
%                 if length(check) >= 70
%                     if leng > 1
%                         check2 = cell2mat(epoch(2));   % GNGSA
%                         check2 = check2(2:6);
%                          if check2 == GSA
%                         switch leng
%                             case 6
%                                 gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                 gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                 gpgsv = epoch(3:6,1);
%                                 nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
%                                 nmeaqm(i,1) = {nmea};
%                             case 5
%                                 gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                 gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                 gpgsv = epoch(3:5,1);
%                                 nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
%                                 nmeaqm(i,1) = {nmea};
%                             case 4
%                                 gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                 gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                 gpgsv = epoch(3:4,1);
%                                 nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
%                                 nmeaqm(i,1) = {nmea};
%                             case 3
%                                 gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                 gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                 gpgsv = epoch(3,1);
%                                 nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
%                                 nmeaqm(i,1) = {nmea};
%                             case 2
%                                 gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                 gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                 gpgsv = [];
%                                 nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
%                                 nmeaqm(i,1) = {nmea};
%                         end
%                     elseif check2 == '$GPGSV'
%                         switch leng
%                             case 5
%                                 gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                 gpgsa = [];
%                                 gpgsv = epoch(2:5,1);
%                                 nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
%                                 nmeaqm(i,1) = {nmea};
%                             case 4
%                                 gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                 gpgsa = [];
%                                 gpgsv = epoch(2:4,1);
%                                 nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
%                                 nmeaqm(i,1) = {nmea};
%                             case 3
%                                 gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                 gpgsa = [];
%                                 gpgsv = epoch(2:3,1);
%                                 nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
%                                 nmeaqm(i,1) = {nmea};
%                             case 2
%                                 gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                 gpgsa = [];
%                                 gpgsv = epoch(2,1);
%                                 nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
%                                 nmeaqm(i,1) = {nmea};
%                         end
%                     end
%                 else
%                     gpgga = strsplit(cell2mat(epoch(1)),{',','*'});
%                     gpgsa = [];
%                     gpgsv = [];
%                     nmea = mknmea(YMD, gpgga, gpgsa, gpgsv);
%                     nmeaqm(i,1) = {nmea};
%                 end
%             end
%             end
% % end
%         case 2
%             for i = 1:length(UBXNMEA(:,1))
%                 epoch = UBXNMEA{i,1};
%                 check = cell2mat(epoch(1));   % GNGGA
%                 leng = length(epoch);
%                 if length(check) >= 70
%                     if leng > 1
%                         check2 = cell2mat(epoch(2));   % GNGSA
%                         check2 = check2(2:6);
%                         if check2 == GSA
%                             gp = cell2mat(epoch(4));   % GPGSV
%                             numgp = str2num(gp(8));
%                             gl = cell2mat(epoch(length(epoch)));    % GLGSV
%                             numgl = str2num(gl(8));
%                             num_line = 3 + numgp + numgl;
%                             switch num_line
%                                 case 11
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                     glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
%                                     gpgsv = epoch(4:3+numgp,1);
%                                     glgsv = epoch(3+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 10
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                     glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
%                                     gpgsv = epoch(4:3+numgp,1);
%                                     glgsv = epoch(3+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 9
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                     glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
%                                     gpgsv = epoch(4:3+numgp,1);
%                                     glgsv = epoch(3+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 8
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                     glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
%                                     gpgsv = epoch(4:3+numgp,1);
%                                     glgsv = epoch(3+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 7
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                     glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
%                                     gpgsv = epoch(4:3+numgp,1);
%                                     glgsv = epoch(3+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 6
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                     glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
%                                     gpgsv = epoch(4:3+numgp,1);
%                                     glgsv = epoch(3+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 5
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                     glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
%                                     gpgsv = epoch(4:3+numgp,1);
%                                     glgsv = epoch(3+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 4
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                     glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
%                                     gpgsv = epoch(4:3+numgp,1);
%                                     glgsv = epoch(3+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 3
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                     glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
%                                     gpgsv = epoch(4:3+numgp,1);
%                                     glgsv = epoch(3+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 2
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                     glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
%                                     gpgsv = epoch(4:3+numgp,1);
%                                     glgsv = epoch(3+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 1
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = strsplit(cell2mat(epoch(2)),{',','*'});
%                                     glgsa = strsplit(cell2mat(epoch(3)),{',','*'});
%                                     gpgsv = epoch(4:3+numgp,1);
%                                     glgsv = epoch(3+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                             end
%                         elseif check2 == 'GPGSV'
%                             switch num_line
%                                 case 11
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = [];
%                                     glgsa = [];
%                                     gpgsv = epoch(2:1+numgp,1);
%                                     glgsv = epoch(1+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 10
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = [];
%                                     glgsa = [];
%                                     gpgsv = epoch(2:1+numgp,1);
%                                     glgsv = epoch(1+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 9
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = [];
%                                     glgsa = [];
%                                     gpgsv = epoch(2:1+numgp,1);
%                                     glgsv = epoch(1+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 8
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = [];
%                                     glgsa = [];
%                                     gpgsv = epoch(2:1+numgp,1);
%                                     glgsv = epoch(1+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 7
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = [];
%                                     glgsa = [];
%                                     gpgsv = epoch(2:1+numgp,1);
%                                     glgsv = epoch(1+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 6
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = [];
%                                     glgsa = [];
%                                     gpgsv = epoch(2:1+numgp,1);
%                                     glgsv = epoch(1+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 5
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = [];
%                                     glgsa = [];
%                                     gpgsv = epoch(2:1+numgp,1);
%                                     glgsv = epoch(1+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                                 case 4
%                                     gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                                     gpgsa = [];
%                                     glgsa = [];
%                                     gpgsv = epoch(2:1+numgp,1);
%                                     glgsv = epoch(1+numgp+1:leng,1);
%                                     nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                                     nmeaqm(i,1) = {nmea};
%                             end
%                         end
%                     else
%                         gga = strsplit(cell2mat(epoch(1)),{',','*'});
%                         gpgsa = [];
%                         glgsa = [];
%                         gpgsv = [];
%                         glgsv = [];
%                         nmea = mkubxnmea(YMD, gga, gpgsa, glgsa, gpgsv, glgsv);
%                         nmeaqm(i,1) = {nmea};
%                     end
%                 end
%             end
% end
%
%
%     nmeaqm = nmeaqm(find(~cellfun(@isempty,nmeaqm)),1);
%
%     i=1;
%     for j = 1:length(nmeaqm)
%         if length(nmeaqm{j}) >3
%             for jj = 1:length(nmeaqm{j}(:,1))
%                 NMEAQM(i,1:length(nmeaqm{j}(1,:))) = nmeaqm{j}(jj,1:length(nmeaqm{j}(1,:)));
%                 i = i + 1;
%             end
%         else
%             NMEAQM(j,1:3) = nmeaqm{j};
%         end
%     end
% end
