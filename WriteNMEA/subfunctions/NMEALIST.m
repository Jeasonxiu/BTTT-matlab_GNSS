function [NMEAlist] = NMEALIST(filename)
%
%
%     clear all
%     filename ='rover_D_160601_6.ubx';
%     filename = 'jprA068a_opG.txt';
% filename = 'jprA214f_s7.txt';
% filename = 'T4_170323_k500.cnb';
% filename = 'DSDC320h.txt';
fid = fopen(filename,'r');
list = textscan(fid,'%s');

%% ubx file에서 NMEA 정보만 추출하는 과정
for i = 1:length(list{1})
%     for i = 2386:2446
line = cell2mat(list{1}(i));
    if length(line) >= 6
        aa = min(findstr(line, '$GP'));
        aaa = min(findstr(line, '$GN'));
        aaaa = min(findstr(line, '$GL'));
        aaaaa = min(findstr(line, '$QZ'));
        aaaaaa = min(findstr(line, '$BD'));
        aaaaaaa = min(findstr(line, '$GB'));
        if ~isempty(aa)
            bb = max(findstr(line, '*'));
            if bb - aa > 6 & bb <= (length(line)-2)
                NMEAlist{i,1} = line(aa:bb+2);
            end
        elseif ~isempty(aaa)
            bb = max(findstr(line, '*'));
            if bb - aaa > 6 & bb <= (length(line)-2)
                NMEAlist{i,1} = line(aaa:bb+2);
            end
        elseif ~isempty(aaaa)
            bb = max(findstr(line, '*'));
            if bb - aaaa > 6 & bb <= (length(line)-2)
                NMEAlist{i,1} = line(aaaa:bb+2);
            end
        elseif ~isempty(aaaaa)
            bb = max(findstr(line, '*'));
            if bb - aaaaa > 6 & bb <= (length(line)-2)
                NMEAlist{i,1} = line(aaaaa:bb+2);
            end
        elseif ~isempty(aaaaaa)
            bb = max(findstr(line, '*'));
            if bb - aaaaaa > 6 & bb <= (length(line)-2)
                NMEAlist{i,1} = line(aaaaaa:bb+2);
            end
        elseif ~isempty(aaaaaaa)
            bb = max(findstr(line, '*'));
            if bb - aaaaaaa > 6 & bb <= (length(line)-2)
                NMEAlist{i,1} = line(aaaaaaa:bb+2);
            end
            
            %                 end
            %             elseif ~isempty(aaa)
            %                 bb = max(findstr(line, '*'));
            %                 if bb - aaa > 6
            %                     NMEAlist{i,1} = line(aa:bb+2);
            %                 end
            %                 NMEAlist{i,1} = line(aa:bb+2);
        end
    end
end

NMEAlist = NMEAlist(find(~cellfun(@isempty,NMEAlist)),1);       % 빈 cell을 없애는 부분
