function N = N_initial(SatsInfo, SatsInfo_before, Obs_Bs, Obs_Rv, prn_RS, N_before)

CCC = 299792458.;  
Freq_L1 = 1575.42e6;    Lambda = CCC/Freq_L1;    % GPS L1��ȣ ���ļ�/����
Sats_compare = SatsInfo - SatsInfo_before;
riselist = find(Sats_compare > 0);
setlist = find(Sats_compare < 0);
N_sd = ((Obs_Bs(:,3) - Obs_Rv(:,3))*Lambda - (Obs_Bs(:,2) - Obs_Rv(:,2))) / Lambda;
% N_sd = ((Obs_Bs(:,3) - Obs_Rv(:,3)) - (Obs_Bs(:,2) - Obs_Rv(:,2))) / Lambda;
%% N_���� ����
if length(riselist) > 0
    for i = 1:length(riselist)
        prn_OS = riselist(i);
        if prn_OS ~= prn_RS
            N_before(prn_OS,1) = N_sd(prn_RS) - N_sd(prn_OS);
        else
            N_before(prn_RS ,1) = 0;
        end
    end
end
%% RS row �� set ���� ����
if length(setlist) > 0
    N_before(setlist) = 0;
    N_before(prn_RS) = 0;
end
N = N_before(find(N_before(:,1) ~= 0),:);