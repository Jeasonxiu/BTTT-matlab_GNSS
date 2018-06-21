function [year month day hour minute second] = utc2gpst(UTC, leap_second)
%
% Do : transform from UTC to GPTS
%
% <input>
%      UTC : Coordinated Universal Time
%            6 by 1 matrix or 1 by 6 matrix format
%      leap_second : 윤초(미입력시:16)
%
% <output>
%      Global Positioning System time
%
% Copyright: taeil Kim, January 7, 2015@INHA University

year    = UTC(1);
month   = UTC(2);
day     = UTC(3);
hour    = UTC(4);
minute  = UTC(5);
second  = UTC(6);

if nargin<2
    leap_second = 16;   % leap seconds GPS is now ahead of UTC by 16 seconds.
                        % now(2015-01-07 15:26:43 +4110)
end

DaysOfMonth = [31,28,31,30,31,30,31,31,30,31,30,31];

if month == 2
    if rem(year,4) == 0
        DaysOfMonth(2) = 29;
    end
end                     % 윤년 확인

[second c] = troops(60, second, leap_second, 0);
[minute c] = troops(60, minute, 0, c);
[hour   c] = troops(24, hour, 0, c);
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