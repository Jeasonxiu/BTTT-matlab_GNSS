function N = N_init(Obs_Bs, Obs_Rv, idx_RS)

CCC = 299792458.;  
Freq_L1 = 1575.42e6;    Lambda = CCC/Freq_L1;    % GPS L1신호 주파수/파장
N_sd = (((Obs_Bs(:,3) - Obs_Rv(:,3))*Lambda - (Obs_Bs(:,2) - Obs_Rv(:,2)))) / Lambda;
%% N_조합 생성
for i = 1:length(Obs_Bs(:,1))
    if i ~= idx_RS
        N_temp(i,1) = N_sd(idx_RS) - N_sd(i);
    else
        N_temp(i,1) = 0;
    end
end
%% RS row 삭제
N = N_temp(find(N_temp(:,1) ~= 0),:);