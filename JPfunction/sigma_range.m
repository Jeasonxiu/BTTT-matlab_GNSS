function sigma = sigma_range(data, a);

% sigma 계산 함수
%
% input data = N X 1 matrix
%       a = sigma 수 (ex> 1 -> 1sigma, 2 -> 2sigma) or % (68.3%, 90%, 95%)
%      
% output two_sigma : 1X2 matrix(2 sigma range)
%
%
% note : 1 sigma(68.3%), 2 sigma (95.5%), 3 sigma (99.7%)
% coded by Joonseong, Feb 6, 2017

STD = std(data);
MEAN = mean(data);

if a < 10
    if a == 1
        sigma = [MEAN-STD, MEAN+STD];
    elseif a == 2
        sigma = [MEAN-2*STD, MEAN+2*STD];
    elseif a == 3
        sigma = [MEAN-3*STD, MEAN+3*STD];
    end
elseif a > 10
    if a == 68
        sigma = [MEAN-STD, MEAN+STD];
    elseif a == 95
        sigma = [MEAN-1.96*STD, MEAN+1.96*STD];
    elseif a == 95.5
        sigma = [MEAN-2*STD, MEAN+2*STD];
    elseif a == 99.7
        sigma = [MEAN-3*STD, MEAN+3*STD];
    end
end
        