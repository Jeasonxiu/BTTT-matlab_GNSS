function [signalOrder, numbers] = GetLineCt(fid)
% function [LeapSec] = GetLeapCt(fid)
%
% 1. Reads one low.
ready = 0;
while ~ready
    s = fgetl(fid);
%     if s(61:68) == 'INTERVAL'
%         interv = str2num(s(5:6));
%     end
    if length(s) > 72
        if s(61:72) == '# / TYPES OF'
            oneLine = s(1:length(s));
            
        end
        
        if s(61:72)=='END OF HEADE'
            ready = 1;
        end
    end
    
end
lineCt=oneLine;
%interval = interv;


% 2. Segmentation of order of the recording.     
% The value of pos{1,1,:} indicates signal character.
% 1 = L1, 2 = L2, 3 = C1, 4 = P1, 5 = P2, 6 = D1, 7 = D2
% The pages of the matrix indicate order of the recording.
% The values of compos is order of each singal charater. zero is without signal.
numbers = sscanf(lineCt,'%d',1);
celltype= {'L1','L2','C1','P1','P2','D1','D2'};
order = 1;
for i=1:7
    str = char(celltype(1,i));
    poscell = findstr(str,lineCt);
    if length(poscell) == 0
       pos{1,1,order} = i;
       pos{1,2,order} = 60; 
    else
       pos{1,1,order} = i;
       pos{1,2,order} = poscell;       
    end    
    order = order + 1;
end
for i=1:6
    for j=i+1:7
    if pos{1,2,i}> pos{1,2,j}
        k1 = pos{1,1,j};
        k2 = pos{1,2,j};
        pos{1,1,j} = pos{1,1,i};
        pos{1,2,j} = pos{1,2,i};
        pos{1,1,i} = k1;
        pos{1,2,i} = k2;
    end
    end
end

compos = zeros(1,7);
for i=1:7
    compos(1,i) = pos{1,1,i};
end
signalOrder = compos;