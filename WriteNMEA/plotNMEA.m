function plotNMEA(NMEAQM, TruePos)

% SYNTAX:
%   plotNMEA(NMEAQM)
%
% INPUT:
%   NMEAQM : after writeNMEA(filename, year, month, day)
%   TurePos : [1X3] coordinates(ECEF)
%
% OUTPUT:
%   Posterror

% test
% clear all
% filename = '362NMEA.txt';
% [NMEAQM] = writeNMEA(filename,2017,12,28);


% TruePos = [-3026676.0349   4067187.8093   3857246.8611];

FinalTTs = unique(NMEAQM(:,1));
NMEA_xyz(:,5:9) = zeros(length(FinalTTs),5);
for i=1:length(FinalTTs)
    gs = FinalTTs(i);
    NMEA_llh = NMEAQM(min(find(NMEAQM(:,1) == gs)), 2:4);
    NumofGPS = length(find(NMEAQM(:,1) == gs & NMEAQM(:,7) < 200));
    if isempty(NumofGPS)
        NumofGPS = 0;end
    NumofBDS = length(find(NMEAQM(:,1) == gs & NMEAQM(:,7) > 200 & NMEAQM(:,7) < 300));
    if isempty(NumofBDS)
        NumofBDS = 0;end
    NumofGLO = length(find(NMEAQM(:,1) == gs & NMEAQM(:,7) > 300 & NMEAQM(:,7) < 400));
    if isempty(NumofGLO)
        NumofGLO = 0;end
    NumofQZSS = length(find(NMEAQM(:,1) == gs & NMEAQM(:,7) > 500));
    if isempty(NumofQZSS)
        NumofQZSS = 0;end
    NumofTotal = length(find(NMEAQM(:,1) == gs));
    if isempty(NumofTotal)
        NumofTotal = 0;end
    NMEA_xyz(i,1) = [gs];
    NMEA_xyz(i,2:4) = gd2xyz(NMEA_llh);
    NMEA_xyz(i,5:9) = [NumofTotal, NumofGPS, NumofBDS, NumofGLO, NumofQZSS];
end


[dXYZ, dNEV] = PosTErrorsNMEA(FinalTTs, TruePos, NMEA_xyz(:,2:4),NMEA_xyz(:,5:9));
