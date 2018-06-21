function [year month day hour minute second] = kst2utc(KST)
%
% Do : transform from KST to UTC
%
% <input>   KST
% <output>  UTC
%
%   Copyright: taeil Kim, February 20, 2015@INHA University

year    = KST(1);
month   = KST(2);
day     = KST(3);
hour    = KST(4);
minute  = KST(5);
second  = KST(6);

DaysOfMonth = [31,31,28,31,30,31,30,31,31,30,31,30];

if month == 2
    if rem(year,4) == 0
        DaysOfMonth(3) = 29;
    end
end                     % 윤년 확인

[hour   c] = troops(24, hour, -9, 0);
[day    c] = troops(DaysOfMonth(month), day-1, 0, c);
[month  c] = troops(12, month-1, 0, c);
[year dum] = troops(9999, year-1, 0, c);

day = day + 1;
month = month + 1;
year = year + 1;

end

%
%
function [output cout] = troops(troop, input, plus, cin)
% troops : (troop)진법 연산
    output=input+plus+cin;
    cout  =floor(output/troop);
    output=mod(output,troop);
end