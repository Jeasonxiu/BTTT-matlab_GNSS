function [GPGSA] = ubxGPGSA(NMEAlist)
    %
    %
    GPGSA = {};
    for i = 1:length(NMEAlist)
        line = NMEAlist{i,1};
        if length(line) >= 6
            check = line(1:6);
            if check == '$GPGSA'
                GPGSA(i,1) = {line};
            else

            end
        end
    end
    if ~isempty(GPGSA)
        GPGSA = GPGSA(find(~cellfun(@isempty,GPGSA)),1);
    end
    