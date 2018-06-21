function [GPSPRC, BDSPRC, GLOPRC, PRC_sorted] = PRCsort(filename, FinalTTs, PRC);
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
% filename = 'JPRT160524.t1';
% QM = load('Q160524_ubx1');
% filename = 'PPS1_170125.t41';
% QM = load('QM170125_A');
if nargin > 2
    if isempty(PRC)
        prc_all = load(filename);
    else
        prc_all = PRC;
    end
end
prc_all(:,1) = prc_all(:,1);
index_gps = find(prc_all(:,2) < 200);
index_bds = find(prc_all(:,2) > 200 & prc_all(:,2) < 300);
index_glo = find(prc_all(:,2) > 300);
gpsprc = prc_all(index_gps,:);
gpsprc(:,2) = gpsprc(:,2);
% gpsprc(:,2) = gpsprc(:,2) - 100;
bdsprc = prc_all(index_bds,:);
bdsprc(:,2) = bdsprc(:,2);
% bdsprc(:,2) = bdsprc(:,2) - 100;
gloprc = prc_all(index_glo,:);
% gloprc(:,2) = gloprc(:,2);
gloprc(:,2) = gloprc(:,2);
% gloprc(:,2) = gloprc(:,2) - 300;
FinalTTs = round(FinalTTs);
if min(gpsprc(:,1)) < min(FinalTTs(:,1))
    gs_min = min(FinalTTs(:,1));
else
    gs_min = min(gpsprc(:,1));
end
if max(gpsprc(:,1)) > max(FinalTTs(:,1))
    gs_max = max(FinalTTs(:,1));
else
    gs_max = max(gpsprc(:,1));
end

all_index_min = find(prc_all(:,1) == gs_min);
all_index_max = find(prc_all(:,1) == gs_max);
PRC_sorted = prc_all(min(all_index_min):max(all_index_max),:);
gps_index_min = find(gpsprc(:,1) == gs_min);
gps_index_max = find(gpsprc(:,1) == gs_max);
bds_index_min = find(bdsprc(:,1) == gs_min);
bds_index_max = find(bdsprc(:,1) == gs_max);
glo_index_min = find(gloprc(:,1) == gs_min);
glo_index_max = find(gloprc(:,1) == gs_max);
GPSPRC = gpsprc(min(gps_index_min):max(gps_index_max),:);
BDSPRC = bdsprc(min(bds_index_min):max(bds_index_max),:);
GLOPRC = gloprc(min(glo_index_min):max(glo_index_max),:);




