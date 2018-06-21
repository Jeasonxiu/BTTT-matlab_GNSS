function x_before = makexbefore(x, prn_RS, SatsInfo_Rv);
k = 1;
for i = 1:32
    prn = i;
    if (SatsInfo_Rv(i,1) == prn) & (SatsInfo_Rv(i,1) ~= prn_RS)
        x_before(i,1) = x(k+3);
        k = k + 1;
    else
        x_before(i,1) = 0;
    end
end

