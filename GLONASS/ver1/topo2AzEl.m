function AzEl = topo2AzEl(topo)

dif_n = topo(1);
dif_e = topo(2);
dif_v = topo(3);

s = sqrt(dif_n^2+dif_e^2+dif_v^2);
AzEl = zeros(1,2);
az = atan2(dif_e,dif_n)*180/pi;

% 방위각은 음수로 나오면 안됨
 if az < 0
     az = az + 360;
end

AzEl(1,1) = az;
AzEl(1,2) = (asin(dif_v/s))*180/pi;