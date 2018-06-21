function [W] = DD_weighting_elsnr(elBs, elRv, SNRBs, SNRRv, W_index)

% SYNTAX:
%   [W] = DD_weighting_elsnr(elBs, elBsv, SNRBs, SNRRv, W_index);
%
% INPUT:
%   elBs = satellite elevations (Base)
%   elBsv = satellite elevations (Rover)
%   SNRBs = signal-to-noise ratio (Base)
%   SNRRv = signal-to-noise ratio (Rover)
%   W_index = coefficient of Weighting 
%
% OUTPUT:
%   W = code-code Weighting
%       parameter used to select the weight mode for GPS observations
%          - W_index=0: same weight for all the observations
%          - W_index=1: weight based on satellite elevation (sin)
%          - W_index=2: weight based on satellite elevation (exp)
%          - W_index=3: weight based on signal-to-noise ratio
%          - W_index=4: weight based on combined elevation and signal-to-noise ratio
%
% Copyright: Joon-seong Gim, PPSoln, April 7, 2017
%

%% Weight Index에 따라서 가중치 계수 설정
%% 0일때, 가중치는 1
snr_a = 30;
snr_0 = 10;
snr_1 = 50;
snr_A = 30;
elea  = 10; % default value for the exponential elevation weight function

if (W_index == 0)

    %code-code or phase-phase co-factor matrix Q construction
    W = 1;

else
    if (W_index == 1)

        %weight vectors (elevation)
        WBs = 1 ./ (sin(elBs * pi/180).^2);
        WRv = 1 ./ (sin(elRv * pi/180).^2);
        
    elseif (W_index == 2)
        %weight vectors (elevation, exponential function)
        eleref = min(elBs)* pi/180; % this is the value for the elevation cut-off angle
        WBs = (1 + elea*exp(-(elBs * pi/180)/eleref)).^2;
        eleref = min(elRv)* pi/180; % this is the value for the elevation cut-off angle
        WRv = (1 + elea*exp(-(elRv * pi/180)/eleref)).^2;

    elseif (W_index == 3)

        %weight vectors (signal-to-noise ratio)
        WBs = 10.^(-(SNRBs-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(SNRBs-snr_1)+1);
        WBs(SNRBs >= snr_1) = 1;
        WRv = 10.^(-(SNRRv-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(SNRRv-snr_1)+1);
        WRv(SNRRv >= snr_1) = 1;

    elseif (W_index == 4)
        %weight vectors (elevation and signal-to-noise ratio)
        WBs = 1 ./ (sin(elBs * pi/180).^2) .* (10.^(-(SNRBs-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(SNRBs-snr_1)+1));
        WBs(SNRBs >= snr_1) = 1;
        WRv = 1 ./ (sin(elRv * pi/180).^2) .* (10.^(-(SNRRv-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(SNRRv-snr_1)+1));
        WRv(SNRRv >= snr_1) = 1;
    end

    %code-code or phase-phase Weighting construction
    W = (WBs + WRv);
end
