function [nmeamsg] = NMEA_extract(NMEAlist,NMEA)
    %
    %
    nmeamsg = {};
    for i = 1:length(NMEAlist)
        line = NMEAlist{i,1};
        if length(line) >= 6
            check = line(1:6);
            if check == NMEA
                nmeamsg(i,1) = {line};
            else

            end
        end
    end
    if ~isempty(nmeamsg)
        nmeamsg = nmeamsg(find(~cellfun(@isempty,nmeamsg)),1);
    end
    