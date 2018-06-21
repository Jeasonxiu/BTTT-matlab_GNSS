function [icol] = PickEPH(eph, sv, time)
%
%function [icol] = PickEPH(eph, sv, time)
%
%   Picks the proper column in the ephemerides array
%               
%       - Originally coded by Per Jarlemark, Harvard-Smithsonian Center for Astrophysics
%       - Modified by Kwan-Dong Park, February, 2002, Harvard-Smithsonian CfA
%

icol = 0;
isat = find(eph(:,18) == sv);

n = length(isat);

if n == 0
    return
end;

icol = isat(1);
dtmin = eph(icol,8) - time;

for k=2:n

    kk = isat(k);
    dt = eph(kk,8) - time;

%    if dt < 0
        if abs(dt) < abs(dtmin)
            icol = kk;
            dtmin = dt;
        end
%    end

end
