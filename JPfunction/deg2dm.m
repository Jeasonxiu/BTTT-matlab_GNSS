function dm = deg2dm(deg)

d = floor(deg);     % deg 
m = mod(deg,1)*60;  % minute

dm = d*100 + m;     % DDMM.MMMMMMMM
