function [dDCB] = AppDCB_glo(DCB, prn)
%
%function [dDCB] = AppDCB_glo(DCB, prn);
%
% DO: Extract and convert GLONASS DCB from the input DCB array 
%
% <input>   DCB: DCB array  -- FORMAT c1(PRN); c2(BIAS[ns]); c3(RMS)
%           prn: GLONASS prn  
%
% <output>  dDCB: DCB in seconds
%
%--- History ---
% 1/4/2015 김미소 버전 AppDCBsat을 살짝 수정함
%

%% 글로나스 위성은 PRN에 200을 더해서 DCB 행렬에 저장되어 있음 
prn = prn + 200;
%% 해당 위성을 찾아서 주파수 고려 및 단위 변환 수행
indx = find(DCB(:,1) == prn);
dDCB = 1.5457*DCB(indx, 2); %: L1/L2에서 L1-only로 변환
dDCB = dDCB*1e-9;           %: ns를 second로 변환