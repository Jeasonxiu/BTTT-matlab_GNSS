function [] = PLOTUBXNMEAQM(NMEAQM);
    %
    %
    lati = NMEAQM(:,2);
    longi = NMEAQM(:,3);
    la = fix(lati./100) + (lati./100-fix(lati./100))*100/60;
    lo = fix(longi./100) + (longi./100-fix(longi./100))*100/60;
    figure(102) % google map plot
    hold on
    plot(lo, la,'yo','markersize',3)

    axis([126.8768 126.8776 37.47905 37.47955]);
    plot_google_map;