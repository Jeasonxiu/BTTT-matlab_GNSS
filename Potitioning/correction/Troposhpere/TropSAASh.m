function [ZHD] = TropSAASh(vec_site)
p = 1013.25;
el = 90;
[h] = xyz2gd(vec_site);                 % height

ZHD = (0.002277 * p * 1000 ) / (1 - 0.00266 * cosd(2 * h(1)) - 0.00028 * h(3) / 1000);