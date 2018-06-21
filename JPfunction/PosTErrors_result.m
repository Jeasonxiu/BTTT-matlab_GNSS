function [dXYZ, dNEV] = PosTErrors_result(ttPlot, TruePos, EstPosT)
%
% function [dXYZ, dNEV] = PosTErrors_result(tt, TruePos, EstPosT)
%
% <input>   ttPlot: Time Tags [gs]
%           TruePos: True Position [1x3]
%           EstPosT: Estimated Positioni and Clock Error; c1-3(XYZ Estimates), c4(Clock Estimates)
%
% <output>  dXYZ: XYZ differences wrt TruePos
%           dNEV: NEV differences wrt TruePos
%
% Copyright: Kwan-Dong Park, 9/20/2014
% Modified by Joonseong, Jan 20, 2016
%
%% 입력받은 자료의 크기 및 참값에 해당하는 위경도 결정 결정
NoPos = length(EstPosT);
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2); 
%% 추정된 XYZ와 참값 XYZ의 차이값 결정
dXYZ = zeros(NoPos,3);
for k = 1:NoPos
    dXYZ(k,:) = EstPosT(k,1:3) - TruePos;
end
%% dXYZ를 dNEV로 변환
dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
%% 각 성분별 RMS 계산
% dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
% dNE = sqrt(dN.^2 + dE.^2);        %rmsH = myRMS(dNE);
%                                   %rmsV = myRMS(dV);
% d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
% fprintf('RMS Values \n H:%8.2f \n V:%8.2f \n 3D:%7.2f\n', rmsH, rmsV, rms3)

