function [GNGGA] = ubxGNGGA(NMEAlist)
    %
    %
    GNGGA = {};
    for i = 1:length(NMEAlist)
        line = NMEAlist{i,1};
        if length(line) >= 6
            check = line(1:6);
            if check == '$GNGGA'
                GNGGA(i,1) = {line};
            else

            end
        end
    end
    if ~isempty(GNGGA)
        GNGGA = GNGGA(find(~cellfun(@isempty,GNGGA)),1);
    end
    