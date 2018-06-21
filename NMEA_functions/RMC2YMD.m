function [YMD] = RMC2YMD(GPRMC)
if length(GPRMC) >= 60
    index = findstr(GPRMC,',');
    yymmdd = GPRMC(index(9)+1:index(10)-1);
    if str2num(yymmdd(5:6)) >= 80
        yyyy = 1900 + str2num(yymmdd(5:6));
    else
        yyyy = 2000 + str2num(yymmdd(5:6));
    end
    mm = str2num(yymmdd(3:4));
    dd = str2num(yymmdd(1:2));
    YMD = [yyyy mm dd];
else
    yyyy = 0; mm = 0; dd = 0;
    YMD = [yyyy mm dd];
end
