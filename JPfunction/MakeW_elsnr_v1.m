function [W] = MakeW_elsnr_v1(elR,snr_R, W_index)
%
%function [W] = MakeW_elsnr_v1(el,S1)
%
% DO: Assign an weight for the given elevation angle [pseudo-range]
%
% <input>   el: elevation angle
%           snr: snr
%           W_index : coefficient of Weighting 
%       parameter used to select the weight mode for GPS observations
%          - W_index=0: same weight for all the observations
%          - W_index=1: weight based on satellite elevation (sin)
%          - W_index=2: weight based on satellite elevation (exp)
%          - W_index=3: weight based on signal-to-noise ratio
%          - W_index=4: weight based on combined elevation and signal-to-noise ratio
%
% <output>  W: Weight 
%
% Copyright: Kwan-Dong Park, Jipyong Space, October 18, 2014
%

%% Weight Index에 따라서 가중치 계수 설정
%% 0일때, 가중치는 1
snr_a = 30;
snr_0 = 10;
snr_1 = 50;
snr_A = 30;
elea  = 10; % default value for the exponential elevation weight function

%total number of visible satellites
n = length(elR);

if (W_index == 0 || (~any(elR) || ~any(snr_R)))
    
    %code-code or phase-phase co-factor matrix Q construction
    W = eye(n);
    
else
    if (W_index == 1)
        
        %weight vectors (elevation)
        q_R = 1 ./ (sin(elR * pi/180).^2);
        
    elseif (W_index == 2)
        %weight vectors (elevation, exponential function)
        eleref = min(elR)* pi/180; % this is the value for the elevation cut-off angle
        q_R = (1 + elea*exp(-(elR * pi/180)/eleref)).^2;
        
    elseif (W_index == 3)
        
        %weight vectors (signal-to-noise ratio)
        q_R = 10.^(-(snr_R-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_R-snr_1)+1);
        q_R(snr_R >= snr_1) = 1;
        
    elseif (W_index == 4)
        %weight vectors (elevation and signal-to-noise ratio)
        q_R = 1 ./ (sin(elR * pi/180).^2) .* (10.^(-(snr_R-snr_1)/snr_a) .* ((snr_A/10.^(-(snr_0-snr_1)/snr_a)-1)./(snr_0-snr_1).*(snr_R-snr_1)+1));
        q_R(snr_R >= snr_1) = 1;
    end
    
    %code-code or phase-phase co-factor matrix Q construction
    W = diag(q_R);
end
