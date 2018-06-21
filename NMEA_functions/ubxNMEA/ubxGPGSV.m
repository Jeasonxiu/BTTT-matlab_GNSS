function [GPGSV] = ubxGPGSV(NMEAlist)
    %
    %
    GPGSV = {};
    for i = 1:length(NMEAlist)
        line = NMEAlist{i,1};
        if length(line) >= 6
            check = line(2:6);
            if check == 'GPGSV'
                GPGSV(i,1) = {line};
            else
                
            end
        end
    end
    i = 0;
    if ~isempty(GPGSV)
        GPGSV = GPGSV(find(~cellfun(@isempty,GPGSV)),1);
        for i = 1:length(GPGSV)
            line = cell2mat(GPGSV(i,1));
            index = findstr(line,',');
            ToNo = str2num(line(index(1)+1:index(2)-1));
            msgNo = str2num(line(index(2)+1:index(3)-1));
            if ToNo == 4 && msgNo == 1
                line2 = cell2mat(GPGSV(i+1,1));
                line3 = cell2mat(GPGSV(i+2,1));
                line4 = cell2mat(GPGSV(i+3,1));
                gpgsv{i,1} = {line; line2; line3; line4};
            elseif ToNo == 3 && msgNo == 1
                line2 = cell2mat(GPGSV(i+1,1));
                line3 = cell2mat(GPGSV(i+2,1));
                gpgsv{i,1} = {line; line2; line3};
            elseif ToNo == 2 && msgNo == 1
                line2 = cell2mat(GPGSV(i+1,1));
                gpgsv{i,1} = {line; line2};
            elseif ToNo == 1 && msgNo == 1
                glgsv{i,1} = {line; line2};
            else
            end
            
        end
    end
    
    GPGSV = gpgsv;
    GPGSV = GPGSV(find(~cellfun(@isempty,GPGSV)),1);
    