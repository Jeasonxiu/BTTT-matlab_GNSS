function t = rtcm_time_corr(pc_time_gs, mz_count, system)
%
%function t = rtcm_time_corr(pc_time_gs, mz_count, system)
%
%   mz_count 를 이용해 정확한 현재 시각을 맞춤
%       <GLONASS에서 윤초 이용 : 현재 윤초(16) 2015.03.06>
%
%   pc_time_gs : 파일에 기록된 현재 pc 시간
%   mz_count   : Modified Z_Count
%   system     : 위성 시스템(eg. 'GLONASS')
%
%   Copyright: taeil Kim, March 6, 2015@INHA University

pc_time_gs = round(pc_time_gs);
mz_count   = round(mz_count*0.6);           % M.Z_COUNT 스케일

if strcmp(system, 'GLONASS')                % GLONASS는 윤초만큼 더해줘야함
    [mz_count dum] = troops(3600, mz_count, 16, 0);
end

t = round( (pc_time_gs - mz_count)./3600 ) * 3600 + mz_count;
if t >= 604800, t = t - 604800; end         % gs max값 초과시
if t < 0, t = t + 604800; end               % gs min값 미만시
end
%
%
function [output cout] = troops(troop, input, plus, cin)
% troops : (troop)진법 연산
    output=input+plus+cin;
    cout  =floor(output/troop);
    output=mod(output,troop);
end