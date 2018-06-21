function [W] = DDMakeW_snrDiff(S1RS,S1OS)
%
%function [W] = MakeW_pr(el,S1)
%
% DO: Assign an weight for the given elevation angle [pseudo-range]
%
% <input>   el: elevation angle
%
% <output>  W: Weight 
%
% Copyright: Kwan-Dong Park, Jipyong Space, October 18, 2014
%

%% 코드-의사거리 가중치 기본은 1미터로 설정
pr_weight = 1.00;

RS_snr_diff = abs(S1RS(1) - S1RS(2)) ;
OS_snr_diff = abs(S1OS(1) - S1OS(2)) ;
RS_snr = -4 * RS_snr_diff;
OS_snr = -4 * OS_snr_diff + 90;
OS_snr_scaleup = 2 * S1OS(2) -30;
if OS_snr_scaleup >= 90
    OS_snr_scaleup = 90;
end

if RS_snr >= 90
    RS_snr = 89;
end
if OS_snr >= 90
    OS_snr = 89;
end

W_RS_snr_diff = cosd(RS_snr);
W_OS_snr_diff = sind(OS_snr)^2;
W_RS_snr = sind(S1RS(2));
% W_OS_snr = sind(S1OS(2));
% W_OS_snr = sind(S1OS(2))^2;
W_OS_snr = sind(OS_snr_scaleup)^3;
% W = W_RS_snr * W_OS_snr * W_OS_el * W_RS_snr_diff * W_OS_snr_diff;
% W = W_OS_snr * W_OS_snr_diff;
% W = W_OS_snr * W_OS_el * W_OS_snr_diff;
% W = W_OS_snr * W_OS_el;
% W = W_OS_snr;
% W = W_OS_snr * W_OS_snr_diff;
W = W_OS_snr_diff;
%% 10/18/2014 김혜인 고도각 가중치
% W = sind(el)^2;