clear all
close all

load('almanac_17255.mat');
load('eph170911.mat');
alm = [alm_gps;alm_bds];

AppPos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
x = [AppPos 1]; x = x';
x_alm = [AppPos 1]; x_alm = x_alm';

gs = 444400;

prns = intersect(alm(:,1), eph(:,18));

for i=1:length(prns)
    prn = prns(i);
    icol = PickEPH(eph, prn, gs);
    STT = GetSTTbrdc(gs, prn, eph, x(1:3)');                % GPS PP
    tc = gs - STT;                                          % GPS PP
    vec_sat = GetSatPosNC(eph, icol, tc);                   % GPS PP
    vec_sat = RotSatPos(vec_sat, STT);                      % GPS PP 지구자전 고려
    SatPos_eph(i,:) = [prn,    vec_sat'];
    STT_alm = GetSTTalm_GC(gs, prn, alm, x_alm(1:3)');
    tc_alm = gs - STT_alm;                                  % GPS alm PP
    vec_sat_alm = GetSatPos_GC_almanac(alm, prn, tc_alm);      % GPS alm PP
    vec_sat_alm = RotSatPos(vec_sat_alm, STT_alm);          % GPS alm PP 지구자전 고려
    SatPos_alm(i,:)=[prn,    vec_sat_alm'];
    
end