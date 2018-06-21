function [TruePos_H] = Apply_H(TruePos,H)
% TruePos = [-3026789.236 4067255.523 3857098.106];
% H = 1.59;
% 
gd = xyz2gd(TruePos);
gd(3) = gd(3)+H;
TruePos_H=gd2xyz(gd);
% 
% TruePos-TruePos_H


% gd = xyz2gd(TruePos);
% % gd(3) = gd(3) + H;
% % TruePos_H=gd2xyz(gd);
% 
% a = xyz2topo([0 0 H], gd(1), gd(2));
% TruePos_H = TruePos + a;

end

%[-3026789.236 4067255.523 3857098.106]
%[-3026789.990 4067256.536 3857099.073]


%[-3026789.236 4067255.523 3857098.106]
%[-3026790.065 4067256.638 3857099.170]