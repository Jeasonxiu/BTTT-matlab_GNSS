function [RMS] = myRMS(value)

% d_vtec=lim(:,2:5)-[gim(:,2) gim(:,2) gim(:,2) gim(:,2)];
% rmse_gr=sqrt(mean(d_vtec(:,1).^2));

RMS = sqrt(mean(value(:,1).^2));