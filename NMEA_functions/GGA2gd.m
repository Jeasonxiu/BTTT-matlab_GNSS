function [gd, gs, utc, la, lo, h, xyz] = GGA2gd(filename)
%
% function [gd, gs, utc, la, lo, h, xyz] = GGA2gd(filename)
%
%   Read the NMEA file, make gd from GPGGA
%
%   input filename : GPGGA after getGPGGA
%
%   Example : [gd, la, lo, h, xyz, utc] = GGA2gd('NMEA.txt')
%
%   coded by Joonseong Gim, Jan 21, 2016
%   GPGGA 내 UTC를 gw, gs로 변환하는 기능 추가

% yy = 2016;, mm = 01; dd = 11;
% filename = 'DAUU049i.ubx';
% filename = '(A)20160111153221.txt';
% filename = 'jprA050i_tab.txt';
% filename = 'hyu0215.ubx';

fid=fopen(filename,'r');
if fid == -1
    disp('Cannot locate the input file!')
    
else
    GPGGA = getGPGGA(filename);
    for i = 1:length(GPGGA)
        line = cell2mat(GPGGA(i,1));
        index = findstr(line,',');
        check = length(line);
        if check >= 70
            epoch{i,1} ={line};
        end
    end
    gpgga = epoch(find(~cellfun(@isempty,epoch)),1);
    
    GPRMC = getGPRMC(filename);
    if ~isempty(GPRMC)
        for ii = 1: length(GPRMC)
            line = cell2mat(GPRMC(ii,1));
            if length(line) >= 70
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
        else
            yyyy = 0;
            mm = 0;
            dd = 0;
        end
    else
        yyyy = 0;
        mm = 0;
        dd = 0;
    end
    
    
    if yyyy > 0
        
        for j = 1: length(gpgga)
            line = cell2mat(gpgga{j}(1));
            index = findstr(line,',');
            utc(j,:) = str2num(line(index(1)+1:index(2)-1));
            utc_h(j,1) = str2num(line(index(1)+1:index(1)+2));
            utc_m(j,1) = str2num(line(index(1)+3:index(1)+4));
            utc_s(j,1) = str2num(line(index(1)+5:index(1)+6));
            la_dms(j,1) = str2num(line(index(2)+1:index(3)-1));
            la(j,1) = fix(la_dms(j,1)/100) + (la_dms(j,1)/100-fix(la_dms(j,1)/100))*100/60;
            lo_dms(j,1) = str2num(line(index(4)+1:index(5)-1));
            lo(j,1) = fix(lo_dms(j,1)/100) + (lo_dms(j,1)/100-fix(lo_dms(j,1)/100))*100/60;
            H(j,1) = str2num(line(index(9)+1:index(10)-1));
            N(j,1) = str2num(line(index(11)+1:index(12)-1));
            h(j,1) = H(j,1) + N(j,1);
        end
        gd = [la lo h];
        
        for iii = 1:length(gd(:,1))
            [gw,Gs] = date2gwgs(yyyy,mm,dd,utc_h(iii),utc_m(iii),utc_s(iii));
            gs(iii,1) = round(Gs);
        end
        %% gd 행렬을 xyz 행렬로 변환하는 과정
        i = 0;
        for i = 1:length(gd)
            xyz(i,:) = gd2xyz(gd(i,:));
            i = i + 1;
        end
    else
        
        for j = 1: length(gpgga)
            line = cell2mat(gpgga{j}(1));
            index = findstr(line,',');
            utc(j,:) = str2num(line(index(1)+1:index(2)-1));
            utc_h(j,1) = str2num(line(index(1)+1:index(1)+2));
            utc_m(j,1) = str2num(line(index(1)+3:index(1)+4));
            utc_s(j,1) = str2num(line(index(1)+5:index(1)+6));
            la_dms(j,1) = str2num(line(index(2)+1:index(3)-1));
            la(j,1) = fix(la_dms(j,1)/100) + (la_dms(j,1)/100-fix(la_dms(j,1)/100))*100/60;
            lo_dms(j,1) = str2num(line(index(4)+1:index(5)-1));
            lo(j,1) = fix(lo_dms(j,1)/100) + (lo_dms(j,1)/100-fix(lo_dms(j,1)/100))*100/60;
            H(j,1) = str2num(line(index(9)+1:index(10)-1));
            N(j,1) = str2num(line(index(11)+1:index(12)-1));
            h(j,1) = H(j,1) + N(j,1);
            gs(j,1) = 0;
        end
        gd = [la lo h];
        
        %% gd 행렬을 xyz 행렬로 변환하는 과정
        i = 0;
        for i = 1:length(gd)
            xyz(i,:) = gd2xyz(gd(i,:));
            i = i + 1;
        end
    end
end
