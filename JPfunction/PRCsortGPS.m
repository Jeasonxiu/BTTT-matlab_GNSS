function [PRC_Sorted] = PRCsortGPS(filename, QM); 
%
% function [PRC_Sorted] = PRCsort(filename, NMEAQM); 
%
% <input>   filename: prc file
%           NMEAQM: NMEAQM file from logged NMEA file
%
% <output>  PRC_Sorted: Sorting PRC using logged NMEA time
%
% Modified by JOON, 26/02/2016
%

% clear all
% filename = 'JPRT160218.t1';
prc_all = load(filename);
prc_all(:,1) = prc_all(:,1);
index_gps = find(prc_all(:,2) < 300);
prc = prc_all(index_gps,:);
prc(:,2) = prc(:,2) - 100;
gs_min = min(QM(:,1));
gs_max = max(QM(:,1));
index_min = find(prc(:,1) == gs_min);
index_max = find(prc(:,1) == gs_max);
PRC_Sorted = prc(min(index_min):max(index_max),:);




    