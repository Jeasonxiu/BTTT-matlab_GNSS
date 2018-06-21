function [] = rmoutl_stat(yrs, dvs, res);

P = polyfit(yrs,dvs,1); 
Rep = std(res); 
% vel_str = num2str(P(1)*10, '%4.1f');
% rep_str = num2str(std(res)*10, '%4.1f');
%% fprintf¿¡ 'mm'Ãß°¡
fprintf('%4.1fmm %4.1fmm \n', P(1)*10, Rep*10)