function [gd, gs, la, lo, h, xyz] = GGARMC2gd(GPGGA, GPRMC)
%
% function [gd, gs, la, lo, h, xyz] = GGARMC2gd(GPGGA, GPRMC)
%
%   Read the GPGGA GPRMC, make gd from GPGGA
%   
%   input GPGGA : GPGGA(cell array) after getGPGGA
%   input GPRMC : GPRMC(cell array) after getGPRMC
%
%   Example : [gd, gs, la, lo, h, xyz] = GGARMC2gd(GPGGA, GPRMC)
%
%   coded by Joonseong Gim, Jan 28, 2016
%   GPGGA 내 UTC와 GPRMC의 DDMMYY를 가지고  gw, gs로 변환하는 기능 추가

% yy = 2016;, mm = 01; dd = 11;
% filename = 'jprA014a.txt';
% filename = '(A)20160111153221.txt';
% GPGGA = getGPGGA(filename);

for i = 1:length(GPGGA)
    line1 = cell2mat(GPGGA(i,1));
    line2 = cell2mat(GPRMC(i,1));
    check = length(line1);
    if check >= 70
        epoch{i,1} ={line1; line2};
    end
end
ggarmc = epoch(find(~cellfun(@isempty,epoch)),1);
for j = 1: length(ggarmc)
%     line = cell2mat(gpgga{j});
    gga = strsplit(cell2mat(ggarmc{j}(1)),{',','*'});
    rmc = strsplit(cell2mat(ggarmc{j}(2)),{',','*'});
    utc = cell2mat(gga(2));  
    utc_h(j,1) = str2num(utc(1:2));
    utc_m(j,1) = str2num(utc(3:4));
    utc_s(j,1) = str2num(utc(5:6));
    la_dms(j,1) = str2double(cell2mat(gga(3)));
    la(j,1) = fix(la_dms(j,1)/100) + (la_dms(j,1)/100-fix(la_dms(j,1)/100))*100/60;
    lo_dms(j,1) = str2double(cell2mat(gga(5)));
    lo(j,1) = fix(lo_dms(j,1)/100) + (lo_dms(j,1)/100-fix(lo_dms(j,1)/100))*100/60;
    H = str2double(cell2mat(gga(10)));
    N = str2double(cell2mat(gga(12)));
    h(j,1) = H + N;
    
    yymmdd = cell2mat(rmc(10));
    if str2num(yymmdd(5:6)) >= 80
        yyyy(j,1) = 1900 + str2num(yymmdd(5:6));
    else
        yyyy(j,1) = 2000 + str2num(yymmdd(5:6));
    end
    mm(j,1) = str2num(yymmdd(3:4));
    dd(j,1) = str2num(yymmdd(1:2));

end

gd = [la lo h];
 
 for i = 1:length(gd(:,1))
     [gw,gs(i,1)] = date2gwgs(yyyy(i),mm(i),dd(i),utc_h(i),utc_m(i),utc_s(i));
 end
 
%% gd 행렬을 xyz 행렬로 변환하는 과정
i = 0;
for i = 1:length(gd)
    xyz(i,:) = gd2xyz(gd(i,:));
    i = i + 1;
end



        
