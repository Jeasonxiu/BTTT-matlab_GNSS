function [new_PRC p_obs] = SatPos(i, eph, r_obs, y, tgps, app_xyz, Ref_Pos, VRS_Pos, RTCM_1,RTCM_2,RTCM_3,RTCM_4,RTCM_5,RTCM_6,RTCM_7,RTCM_8,RTCM_9,RTCM_10,RTCM_11,RTCM_12,RTCM_13,RTCM_14,RTCM_15)
% function [PRC_old PRC p_obs RMS] = SatPos(eph, r_obs, y, tgps, app_xyz, Ref_Pos, VRS_Pos, RTCM_1,RTCM_2,RTCM_3,RTCM_4,RTCM_5,RTCM_6)
% function [p_obs] = SatPos(eph, r_obs, y, tgps, app_xyz, Ref_Pos, VRS_Pos, RTCM);

% Usage: [p_obs] = SatPos(eph, r1_obs1, prn_vec, tgps, ref_pos)
% 위성의 위치계산 및 신호전달시간 고려
% [Hye-In Kim, Inha University, May 2009.]

% 빛의 속도 
c = 299792458.;

% PRC 산출
for m=2:size(y,1)
    [PRC_old, new_PRC, RMS, con_num] = extract_PRC(y(m,1),tgps,Ref_Pos,VRS_Pos,RTCM_1,RTCM_2,RTCM_3,RTCM_4,RTCM_5,RTCM_6,RTCM_7,RTCM_8,RTCM_9,RTCM_10,RTCM_11,RTCM_12,RTCM_13,RTCM_14);
%     [PRC_old, new_PRC, RMS, con_num] = extract_PRC(y(m,1),tgps,Ref_Pos,VRS_Pos,RTCM_1,RTCM_2,RTCM_3,RTCM_4,RTCM_5,RTCM_6,RTCM_7,RTCM_8,RTCM_9,RTCM_10,RTCM_11,RTCM_12,RTCM_13,RTCM_14,RTCM_15,RTCM_16,RTCM_17);
%     [PRC_old, new_PRC, RMS, con_num] = extract_PRC(y(m,1),tgps,Ref_Pos,VRS_Pos,RTCM_1,RTCM_2,RTCM_3);
%     [PRC_old, PRC, RMS] = ComPRC(y(m,1),tgps,Ref_Pos,VRS_Pos,RTCM_1,RTCM_2,RTCM_3,RTCM_4,RTCM_5,RTCM_6);
%     PRC = PickPRC(RTCM, y(m,1), tgps);
    r_obs(m-1,8) = new_PRC;
end

t_delay = zeros(length(r_obs),1);
count = 0;

% t_gps 시간에 맞는 위성 위치(신호전달시간 고려) 및 PRC 계산
while count < 5
    count = count + 1;
    for k = 2:size(y,1)
        eph_i = PickEPH(eph, y(k,1), tgps);
        if eph_i ~= 0;
            sat_Pos = GetSatPos(eph, eph_i, tgps, t_delay(k-1));
            diffPos = sat_Pos(1:3) - app_xyz;
            distPR(k-1) = norm(diffPos);
            r_obs(k-1,2:7) = [eph(eph_i,18) sat_Pos(1) sat_Pos(2) sat_Pos(3) sat_Pos(4) tgps];
        else
            sat_Pos = [0 0 0 0];
            r_obs(k-1,2:7) = [0 sat_Pos(1) sat_Pos(2) sat_Pos(3) sat_Pos(4) tgps];
        end
    end
    t_delay_new = distPR / c;
    t_delay = t_delay_new;
    p_obs = Remzero(r_obs);
end
