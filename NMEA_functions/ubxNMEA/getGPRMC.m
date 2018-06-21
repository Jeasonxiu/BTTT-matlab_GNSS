function [RMC] = getGPRMC(filename);

% clear all
% filename = 'COM15_161005_032031.ubx';
[NMEAlist] = NMEALIST(filename);
[GNRMC] = ubxGNRMC(NMEAlist);
cnt=0;
for i = 1:length(GNRMC)
    line = GNRMC{i};
    rmc = textscan(line(1:end-4),'%s %s %s %f %s %f %s %f %f %s %f %f','delimiter',',');
    if cell2mat(rmc(4)) ~= 0
        cnt = cnt+1;
        lati = cell2mat(rmc(4));
        degLat = fix(lati/100);
        minLat = lati-degLat*100;
        la = degLat + minLat/60;        % Latitude
        longi = cell2mat(rmc(6));
        degLon = fix(longi/100);
        minLon = longi-degLon*100;
        lo = degLon + minLon/60;        % Latitude
        UTC = rmc(2);
        UTC=cell2mat(UTC{1});
        hh = str2double(UTC(1:2));
        mm = str2double(UTC(3:4));
        ss = str2double(UTC(5:end));
        date = cell2mat(rmc{10});
        day = str2double(date(1:2));
        month = str2double(date(3:4));
        year = str2double(date(5:6));
        if year > 80
            year = year + 1900;
        else
            year = year + 2000;
        end
        [gw, gs] = date2gwgs(year, month, day, hh, mm, ss);
    end
    RMC(cnt,1:3) = [round(gs), la, lo];
        
end