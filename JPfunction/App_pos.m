function [apppos] = App_pos(obs_file)

%2006-05 김상호
%헤더에서 Approx Position 정보 읽어오기

fid = fopen(obs_file,'r');

while 1
    s = fgets(fid);
    if length(s) > 71
        if s(61:72) == 'APPROX POSIT'
            apppos = [str2num(s(1:14)) str2num(s(16:29)) str2num(s(30:43))];
            break
        end
    end
end