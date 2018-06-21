function [GPRMC] = ubxGPRMC(NMEAlist)
    %
    %
    GPRMC = {};
    for i = 1:length(NMEAlist)
        line = NMEAlist{i,1};
        if length(line) >= 6
            check = line(1:6);
            if check == '$GPRMC'
                GPRMC(i,1) = {line};
            else

            end
        end
    end
    if ~isempty(GPRMC)
        GPRMC = GPRMC(find(~cellfun(@isempty,GPRMC)),1);
    end
    