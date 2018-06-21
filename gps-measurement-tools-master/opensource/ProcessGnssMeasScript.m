clear all;
close all;

%ProcessGnssMeasScript.m, script to read GnssLogger output, compute and plot:
% pseudoranges, C/No, and weighted least squares PVT solution
%
% you can run the data in pseudoranges log files provided for you: 
% prFileName = 'pseudoranges_log_2016_06_30_21_26_07.txt'; %with duty cycling, no carrier phase
% prFileName = 'pseudoranges_log_2016_08_22_14_45_50.txt'; %no duty cycling, with carrier phase
% prFileName = 'pseudoranges_log_2017_03_02_15_33_59.txt'; %no duty cycling, with carrier phase
% prFileName = 'pseudoranges_log_2017_08_01_15_39_28.txt';        % 대성 B point 수평
% prFileName = 'pseudoranges_log_2017_08_01_15_50_15.txt';        % 대성 B point 뒷면
% prFileName = 'pseudoranges_log_2017_08_01_16_02_35.txt';        % 대성 B point 세로 수직
% prFileName = 'pseudoranges_log_2017_08_01_16_14_00.txt';        % 대성 B point 가로 수직
% prFileName = 'pseudoranges_log_2017_09_14_16_31_57.txt';      % 대성 B 지점 1
% prFileName = 'pseudoranges_log_2017_09_14_16_42_37.txt';      % 대성 B 지점 2
% prFileName = 'pseudoranges_log_2017_09_14_16_53_37.txt';      % 대성 B 지점 3
% prFileName = 'pseudoranges_log_2017_09_14_17_04_08.txt';      % 대성 B 지점 4
% prFileName = 'pseudoranges_log_2017_09_14_17_14_31.txt';      % 대성 B 지점 5
% prFileName = 'pseudoranges_log_2017_09_14_17_25_10.txt';      % 대성 B 지점 6
% prFileName = 'pseudoranges_log_2017_09_14_17_36_04.txt';      % 대성 B 지점 7
% prFileName = 'pseudoranges_log_2017_09_14_17_46_35.txt';      % 대성 B 지점 8
% prFileName = 'pseudoranges_log_2017_09_14_17_57_01.txt';      % 대성 B 지점 9
% prFileName = 'pseudoranges_log_2017_09_14_18_08_57.txt';      % 대성 B 지점 10
% prFileName = 'pseudoranges_log_2017_10_13_14_46_27.txt';
prFileName = 'S8_PA_gnss_log_2018_05_15_15_14_46.txt';
prFileName = 'gnss_log_2018_05_27_17_10_57.txt';
% as follows
% 1) copy everything from GitHub google/gps-measurement-tools/ to 
%    a local directory on your machine
% 2) change 'dirName = ...' to match the local directory you are using:
% dirName = '~/Documents/MATLAB/gpstools/opensource/demoFiles';
% dirName = 'E:\demoFiles';
% dirName = 'E:\PPSoln\PPSoln-js\matlab\codes\gps-measurement-tools-master\opensource\demoFiles';
% dirName = 'E:\PPSoln\PPSoln-js\matlab\data\Smartphone\170914';
% dirName = 'E:\PPSoln\PPSoln-js\matlab\data\Smartphone\171013';
dirName = 'E:\PPSoln\PPSoln-js\matlab\data\testdata\RTK\180515_JAVAD'
dirName = 'E:\PPSoln\PPSoln-js\matlab\data\testdata\RTK\180527_JAVAD_S8';
% dirName = '/demoFiles/';
% 3) run ProcessGnssMeasScript.m script file 
param.llaTrueDegDegM = [];

%Author: Frank van Diggelen
%Open Source code for processing Android GNSS Measurements

%% data
%To add your own data:
% save data from GnssLogger App, and edit dirName and prFileName appropriately
%dirName = 'put the full path for your directory here';
%prFileName = 'put the pseuoranges log file name here';

%% parameters
%param.llaTrueDegDegM = [];
%enter true WGS84 lla, if you know it:
% param.llaTrueDegDegM = [37.422578, -122.081678, -28];%Charleston Park Test Site
% param.llaTrueDegDegM = [37.4797970132665, 126.876946215514, 157.874877604656];%PPsoln A point Coordinates
param.llaTrueDegDegM = [37.4797132062492,126.876985221076,157.888893392868];    % 대성 B point
%% Set the data filter and Read log file
dataFilter = SetDataFilter;
[gnssRaw,gnssAnalysis] = ReadGnssLogger(dirName,prFileName,dataFilter);
if isempty(gnssRaw), return, end

%% Get online ephemeris from Nasa ftp, first compute UTC Time from gnssRaw:
fctSeconds = 1e-3*double(gnssRaw.allRxMillis(end));
utcTime = Gps2Utc([],fctSeconds);
allGpsEph = GetNasaHourlyEphemeris(utcTime,dirName);
if isempty(allGpsEph), return, end

%% process raw measurements, compute pseudoranges:
[gnssMeas] = ProcessGnssMeas(gnssRaw);

%% plot pseudoranges and pseudorange rates
h1 = figure;
[colors] = PlotPseudoranges(gnssMeas,prFileName);
h2 = figure;
PlotPseudorangeRates(gnssMeas,prFileName,colors);
h3 = figure;
PlotCno(gnssMeas,prFileName,colors);

%% compute WLS position and velocity
% gpsPvt = GpsWlsPvt(gnssMeas,allGpsEph);
% [gpsPvt]= GpsWlsPvt(gnssMeas,allGpsEph);
[gpsPvt,prs,xo, xHat, wpr,svXyzTrx,ttx]= GpsWlsPvt(gnssMeas,allGpsEph);

%% plot Pvt results
h4 = figure;
ts = 'Raw Pseudoranges, Weighted Least Squares solution';
PlotPvt(gpsPvt,prFileName,param.llaTrueDegDegM,ts); drawnow;
h5 = figure;
PlotPvtStates(gpsPvt,prFileName);

%% Plot Accumulated Delta Range 
if any(any(isfinite(gnssMeas.AdrM) & gnssMeas.AdrM~=0))
    [gnssMeas]= ProcessAdr(gnssMeas);
    h6 = figure;
    PlotAdr(gnssMeas,prFileName,colors);
    [adrResid]= GpsAdrResiduals(gnssMeas,allGpsEph,param.llaTrueDegDegM);drawnow
    h7 = figure;
    PlotAdrResids(adrResid,gnssMeas,prFileName,colors);
end
%% end of ProcessGnssMeasScript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2016 Google Inc.
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
