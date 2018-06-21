function [PRC] = PickPRC(rtcm, sv, tgps)
% rtcm = RTCM2; sv = 210; tgps = 126400;
% time=mod(tgps,86400);
% tgps = tgps+14;
% icol = rtcm(:,3) == sv & rtcm(:,1)==tgps;
% PRC = rtcm(icol,4);
% if PRC <= 0
%     PRC = PRC;
% else
%     PRC=0;
% end
tgps=round(tgps);
icol = rtcm(:,2) == sv & rtcm(:,1)==tgps;
PRC = rtcm(icol,3);
if PRC <= 0
    PRC = PRC;
else
    PRC=0;
end