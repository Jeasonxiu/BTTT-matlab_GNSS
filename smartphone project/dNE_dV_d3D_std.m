load('.mat');


[dXYZ, dNEV] = PosTErrors4(estm(:,1), TruePos, estm(:,2:5),visiSat);
[dXYZ_d, dNEV_d] = PosTErrors4(estm_d(:,1), TruePos, estm_d(:,2:5),visiSat);

save(strcat(FileQM,'.mat'), estm, estm_d, TruePos, visiSat)
%% Standard deviation
dNEstd= std(sqrt(dNEV(:,1).^2+dNEV(:,2).^2))
dVstd= std(dNEV(:,3))
dNEVstd= std(sqrt(dNEV(:,1).^2+dNEV(:,2).^2+dNEV(:,3).^2))
dNEstd_d= std(sqrt(dNEV_d(:,1).^2+dNEV_d(:,2).^2))
dVstd_d= std(dNEV_d(:,3))
dNEVstd_d= std(sqrt(dNEV_d(:,1).^2+dNEV_d(:,2).^2+dNEV_d(:,3).^2))