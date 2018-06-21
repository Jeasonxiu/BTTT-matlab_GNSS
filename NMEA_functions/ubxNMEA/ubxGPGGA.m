function [GPGGA] = ubxGPGGA(NMEAlist)
    
    %
    GPGGA = {};
    for i = 1:length(NMEAlist)
        line = NMEAlist{i,1};
        if length(line) >= 6
            check = line(1:6);
            if check == '$GPGGA'
                GPGGA(i,1) = {line};
            else

            end
        end
    end
    if ~isempty(GPGGA)
        GPGGA = GPGGA(find(~cellfun(@isempty,GPGGA)),1);
    end
    