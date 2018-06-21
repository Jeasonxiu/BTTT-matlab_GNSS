% Compute Val = || Xs - X || + b and its Jacobian.
function [Val, Jacob] = PseudorangeEquation(X, SV)

% Each row of SV is the coordinate of a satellite.
% dX = bsxfun(@minus, X([1,3,5])', SV);% X - Xs % 원본
% Val = sum(dX .^2, 2) .^0.5 + X(7);        % 원본
% 
% 
% Jacob = zeros(size(SV, 1), size(X, 1));     % 원본
% Jacob(:, [1,3,5]) = bsxfun(@rdivide, dX, Val);
% Jacob(:, 7) = 1;            % 원본


dX = bsxfun(@minus, X([1,3,5])', SV(:,1:3));% X - Xs
for i=1:length(dX(:,1))
    prn = SV(i,4);
    if prn > 200
        Val(i,:) = sum(dX(i,1:3).^2, 2) .^0.5 + X(8);
    else
        Val(i,:) = sum(dX(i,1:3).^2, 2) .^0.5 + X(7);
    end
end
Jacob = zeros(size(SV, 1), size(X, 1));     % 원본
Jacob(:, [1,3,5]) = bsxfun(@rdivide, dX, Val);
Jacob(:, 7) = 0;
Jacob(:, 8) = 1;
Jacob(:, 9) = 0;
Jacob(:, 10) = 1;
end