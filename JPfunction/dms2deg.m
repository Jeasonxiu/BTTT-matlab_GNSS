function [deg] = dms2deg(dms)


deg = dms(1);
min = dms(2);
sec = dms(3);

deg = deg + min/60 + sec/3600;
fprintf('%13.8f\n',deg);
