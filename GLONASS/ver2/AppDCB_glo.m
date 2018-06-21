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
% 1/4/2015 ��̼� ���� AppDCBsat�� ��¦ ������
%

%% �۷γ��� ������ PRN�� 200�� ���ؼ� DCB ��Ŀ� ����Ǿ� ���� 
prn = prn + 200;
%% �ش� ������ ã�Ƽ� ���ļ� ��� �� ���� ��ȯ ����
indx = find(DCB(:,1) == prn);
dDCB = 1.5457*DCB(indx, 2); %: L1/L2���� L1-only�� ��ȯ
dDCB = dDCB*1e-9;           %: ns�� second�� ��ȯ