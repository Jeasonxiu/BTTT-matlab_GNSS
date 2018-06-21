function t = rtcm_time_corr(pc_time_gs, mz_count, system)
%
%function t = rtcm_time_corr(pc_time_gs, mz_count, system)
%
%   mz_count �� �̿��� ��Ȯ�� ���� �ð��� ����
%       <GLONASS���� ���� �̿� : ���� ����(16) 2015.03.06>
%
%   pc_time_gs : ���Ͽ� ��ϵ� ���� pc �ð�
%   mz_count   : Modified Z_Count
%   system     : ���� �ý���(eg. 'GLONASS')
%
%   Copyright: taeil Kim, March 6, 2015@INHA University

pc_time_gs = round(pc_time_gs);
mz_count   = round(mz_count*0.6);           % M.Z_COUNT ������

if strcmp(system, 'GLONASS')                % GLONASS�� ���ʸ�ŭ ���������
    [mz_count dum] = troops(3600, mz_count, 16, 0);
end

t = round( (pc_time_gs - mz_count)./3600 ) * 3600 + mz_count;
if t >= 604800, t = t - 604800; end         % gs max�� �ʰ���
if t < 0, t = t + 604800; end               % gs min�� �̸���
end
%
%
function [output cout] = troops(troop, input, plus, cin)
% troops : (troop)���� ����
    output=input+plus+cin;
    cout  =floor(output/troop);
    output=mod(output,troop);
end