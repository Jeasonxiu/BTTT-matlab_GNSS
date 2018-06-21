function [GLGSV] = ubxGLGSV(NMEAlist)
    %
    %
    GLGSV = {};
    for i = 1:length(NMEAlist)
        line = NMEAlist{i,1};
        if length(line) >= 6
            check = line(1:6);
            if check == '$GLGSV'
                GLGSV(i,1) = {line};
            else
                
            end
        end
    end
    i=0;
    if ~isempty(GLGSV)
        GLGSV = GLGSV(find(~cellfun(@isempty,GLGSV)),1);
        for i = 1:length(GLGSV)
            line = cell2mat(GLGSV(i,1));
            index = findstr(line,',');
            ToNo = str2num(line(index(1)+1:index(2)-1));
            msgNo = str2num(line(index(2)+1:index(3)-1));
            if ToNo == 4 && msgNo == 1
                line2 = cell2mat(GLGSV(i+1,1));
                line3 = cell2mat(GLGSV(i+2,1));
                line4 = cell2mat(GLGSV(i+3,1));
                glgsv{i,1} = {line; line2; line3; line4};
            elseif ToNo == 3 && msgNo == 1
                line2 = cell2mat(GLGSV(i+1,1));
                line3 = cell2mat(GLGSV(i+2,1));
                glgsv{i,1} = {line; line2; line3};
            elseif ToNo == 2 && msgNo == 1
                line2 = cell2mat(GLGSV(i+1,1));
                glgsv{i,1} = {line; line2};
            elseif ToNo == 1 && msgNo == 1
                glgsv{i,1} = {line};
            else
            end
            
        end
    end
    
    GLGSV = glgsv;
    GLGSV = GLGSV(find(~cellfun(@isempty,GLGSV)),1);
    