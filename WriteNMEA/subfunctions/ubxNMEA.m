function [UBXNMEA] = ubxNMEA(NMEAlist);

for i = 1:length(NMEAlist)
    line = NMEAlist{i};
    if line(2:6) == 'GPGGA' | line(2:6) == 'GNGGA'
        ggaidx(i,1) = i;
    elseif line(2:6) == 'GNRMC' | line(2:6) == 'GPRMC'
        ggaidx(i,1) = i;
    end
end
ggaidx = ggaidx(find(ggaidx(:,1) > 0));
for i = 1:length(ggaidx)-2
    nmea{i,1} = NMEAlist(ggaidx(i):ggaidx(i+1)-1);
    for j = 1:length(nmea{i})
        line = cell2mat(nmea{i}(j));
        if line(2:6) == 'GPGGA' | line(2:6) == 'GNGGA' | line(2:6) == 'GPGSA' | line(2:6) == 'GNGGA'...
                | line(2:6) == 'GPGSV' | line(2:6) == 'GLGSV' | line(2:6) == 'BDGSV' | line(2:6) == 'GPRMC' | line(2:6) == 'GNRMC'...
                | line(2:6) == 'QZGSA'
            NMEA{i,1}(j,1) = {line};
        end
    end
    NMEA{i,1} = NMEA{i,1}(find(~cellfun(@isempty,NMEA{i,1})),1);
end
UBXNMEA = NMEA;