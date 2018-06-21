function [Alpha, Beta] = Get_AB(fid)

ready = 0;

Alpha = zeros(4,1);
Beta = zeros(4,1);

while ~ready
    s = fgets(fid);
        if s(61:67) == 'ION ALP'
            Alpha(1) = str2num(s(4:14));
            Alpha(2) = str2num(s(16:26));
            Alpha(3) = str2num(s(28:38));
            Alpha(4) = str2num(s(40:50));
        end
        if s(61:67) == 'ION BET'
            Beta(1) = str2num(s(1:14));
            Beta(2) = str2num(s(17:26));
            Beta(3) = str2num(s(28:38));
            Beta(4) = str2num(s(40:50));
        end
        if s(61:67)=='END OF '
            ready = 1;
        end
end
