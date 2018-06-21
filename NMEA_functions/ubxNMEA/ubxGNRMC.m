function [GNRMC] = ubxGNRMC(NMEAlist)
    %
    %
    GNRMC = {};
    cnt=1;
    for i = 1:length(NMEAlist)
        line = NMEAlist{i,1};
        if length(line) >= 6
            check = line(1:6);
            if check == '$GNRMC' | check == '$GPRMC'
                GNRMC(cnt,1) = {line};
                cnt=cnt+1;
            else
                
            end
        end
    end
    if ~isempty(GNRMC)
        GNRMC = GNRMC(find(~cellfun(@isempty,GNRMC)),1);
    end
    