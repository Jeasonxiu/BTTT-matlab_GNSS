function rtklib = rtklib_result_conv(filename);
%
%
% input filename : txt file, after RTKlib estimation
% output rtklib : N X 5 matrix(gs, lati, longi, height, quality]
%
% filename ='rtklib_rover1.pos';

fid=fopen(filename,'r');
ready = 0;
while ready == 0
    s = fgets(fid);
    if length(s) > 141
        if s(136:140) == 'ratio'
            ready = 1;
        end
    end
end
ready = 0;
rtklib_rv1 = textscan(fid,'%s'); cnt=1;

for i=1:length(rtklib_rv1{1})
    line = cell2mat(rtklib_rv1{1}(i));
    if length(line) >= 7
        if line(5) =='/'
            year = str2num(line(1:4));
            mo = str2num(line(6:7));
            dd = str2num(line(9:10));
        elseif line(3) == ':'
            hour = str2num(line(1:2));
            mm = str2num(line(4:5));
            sec = str2num(line(7:10));
            [gw gs] = date2gwgs(year, mo, dd, hour, mm, sec);
            rtklib(cnt,1) = [round(gs)];
        elseif line(3) == '.'
            if length(line) == 12
                lati = str2num(line);
                longi = str2num(cell2mat(rtklib_rv1{1}(i+1)));
                height = str2num(cell2mat(rtklib_rv1{1}(i+2)));
                quality = str2num(cell2mat(rtklib_rv1{1}(i+3)));
                ns = str2num(cell2mat(rtklib_rv1{1}(i+4)));
                rtklib(cnt,2:6) = [lati, longi, height, quality, ns];
                cnt = cnt + 1;
            end
        end
    end
end
