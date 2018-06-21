function [nmeaqm] = mknmea2(YMD, gpgga, gngsa, gpgsv, glgsv)
%
% function [nmeaqm] = mknmea(YMD, gpgga, gngsa, gpgsv)
%
%   Read the gpgga, gpgsv, gpgsv per epoch, make a nmea matrix from input
%   
%   input gpgga : gpgga/epoch
%   input gprmc : gprmc/epoch
%   input gngsa : gngsa/epoch
%   input gpgsv : gpgsv/epoch
%
%   output nmeaqm : Num of SA by 8 matrix(gs, lati, longi, alti, prn, el,
%   az, snr)
%
%   Example : [nmeaqm] = mknmea(gpgga, gprmc, gngsa, gpgsv)
%
%   coded by Joonseong Gim, Jan 26, 2016


%% GPRMC 유무에 따른 UTC, gs 선택 과정
%% GGA와 RMC의 UTC, 날짜 정보를 가지고 gw, gs를 만드는 과정
if YMD(1) > 0  % yymmdd 와 UTC 정보를 가지고 gw, gs를 만드는 과정
    utc = cell2mat(gpgga(2)); 
    utc_h = str2num(utc(1:2));
    utc_m = str2num(utc(3:4));
    utc_s = str2num(utc(5:6));
    yyyy = YMD(1); mm = YMD(2); dd =YMD(3);
    [gw,Gs] = date2gwgs(yyyy,mm,dd,utc_h,utc_m,utc_s);
    gs = round(Gs)+17;

    if ~isempty(gngsa)
        gngsa1 = strsplit(cell2mat(gngsa(1)),{',','*'}); gngsa2 = strsplit(cell2mat(gngsa(2)),{',','*'});
        gngsa = [gngsa1(4:length(gngsa1)-4), gngsa2(4:length(gngsa2)-4)];
        if ~isempty(glgsv)

            nmeaqm = zeros(length(gngsa),7);
            nmeaqm(:,1) = gs;
            nmeaqm(:,2) = str2num(cell2mat(gpgga(3)));      % Latitude
            nmeaqm(:,3) = str2num(cell2mat(gpgga(5)));      % Longitude
            nmeaqm(:,4) = str2num(cell2mat(gpgga(10))) + str2num(cell2mat(gpgga(12)));      % Altitude
%             nmeaqm(:,4) = str2num(cell2mat(gpgga(10)));      % Altitude
            GPGSV = GSVmat(gpgsv); GLGSV = GSVmat(glgsv);
            GSV = [GPGSV;GLGSV];
            
            for i = 1:length(gngsa)
                nmeaqm(i,5) = str2num(cell2mat(gngsa(i)));
            end
            [a b] = size(GSV);
            i = 1;
            for j = 1 : a
                for k = 4:4:length(GSV(j,:))
                    SVs(i,1:4) = GSV(j,k:k+3);
                    i = i + 1;
                end
            end
            for ii = 1 : length(nmeaqm(:,5))
                sa = nmeaqm(ii,5);
                for jj = 1 : length(SVs)
                    sv = SVs(jj,1);
                    if sa == sv
                        nmeaqm(ii,6:8) = SVs(jj,2:4);
                    end
                end
            end
        else
            nmeaqm = zeros(length(gngsa),4);
            nmeaqm(:,1) = gs;
            nmeaqm(:,2) = str2num(cell2mat(gpgga(3)));      % Latitude
            nmeaqm(:,3) = str2num(cell2mat(gpgga(5)));      % Longitude
            nmeaqm(:,4) = str2num(cell2mat(gpgga(10))) + str2num(cell2mat(gpgga(12)));      % Altitude
%             nmeaqm(:,4) = str2num(cell2mat(gpgga(10)));      % Altitude
            GPGSV = GSVmat(gpgsv);
            for i = 1:length(gngsa(4:length(gngsa)-4))
                nmeaqm(i,5) = str2num(cell2mat(gngsa(i+3)));
            end
        end
    else
         nmeaqm = zeros(length(gngsa(4:length(gngsa)-4)),3);
         nmeaqm(:,1) = gs;
         nmeaqm(:,2) = str2num(cell2mat(gpgga(3)));      % Latitude
         nmeaqm(:,3) = str2num(cell2mat(gpgga(5)));      % Longitude
         nmeaqm(:,4) = str2num(cell2mat(gpgga(10))) + str2num(cell2mat(gpgga(12)));      % Altitude
%          nmeaqm(:,4) = str2num(cell2mat(gpgga(10)));      % Altitude
    end
else
    utc = cell2mat(gpgga(2)); 
    nmeaqm(:,1) = str2num(utc);
    nmeaqm(:,2) = str2num(cell2mat(gpgga(3)));      % Latitude
    nmeaqm(:,3) = str2num(cell2mat(gpgga(5)));      % Longitude
    nmeaqm(:,4) = str2num(cell2mat(gpgga(10))) + str2num(cell2mat(gpgga(12)));      % Altitude
%     nmeaqm(:,4) = str2num(cell2mat(gpgga(10)));      % Altitude
end


    

