function [GNGSA] = ubxGNGSA(NMEAlist)
%
%
GNGSA = {};
for i = 1:length(NMEAlist)
    line = NMEAlist{i,1};
    if length(line) >= 6
        check = line(1:6);
        if check == '$GNGSA'
            GNGSA(i,1) = {line};
        else
%             GNGSA = {};
        end
    end
end
i = 0;
if ~isempty(GNGSA)
    GNGSA = GNGSA(find(~cellfun(@isempty,GNGSA)),1);
    for i = 1:2:length(GNGSA)
        line = cell2mat(GNGSA(i,1));
        line2 = cell2mat(GNGSA(i+1,1));
        gngsa{i,1} = {line; line2; 0};
    end
    %         for i = 1:length(GNGSA)
    %             line = cell2mat(GNGSA(i,1));
    %             prn = line(12:13);
    %             if prn == ',,'
    %                 line2 = cell2mat(GNGSA(i+1,1));
    %                 gngsa{i,1} = {line; line2; 0};
    %             else
    %                 prn = str2num(prn);
    %                 if prn < 33
    %                     line2 = cell2mat(GNGSA(i+1,1));
    %                     gngsa{i,1} = {line; line2};
    %
    %                 else
    %                 end
    %             end
    %         end
end
if ~isempty(GNGSA)
    GNGSA = gngsa;
    GNGSA = GNGSA(find(~cellfun(@isempty,GNGSA)),1);
end

