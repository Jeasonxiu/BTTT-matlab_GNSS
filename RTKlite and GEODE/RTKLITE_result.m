clear all;
close all;
clc

tic

%% RTKLITE logged file
% filename = 'nmea_20161113.txt';
% filename = 'rtklite_161117.txt';
% filename = 'GGA_2_20170102173904_20170103105410.txt';
% filename = 'GGA_2_20170201.txt';
filename = 'GGA_2_20170124.txt';

%% 실험 날짜 입력
yy = 2017; mo = 1; dd = 24;
%% Reference Point
% TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
% TruePos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
% TruePos = [-3027413.18438068 4047615.83589245 3877061.26661548];    % 20161113
% TruePos = [-3027409.89954179 4047618.47439957 3877061.37645965];    % 20161117
% TruePos = [-3026674.48638231 4067187.35368425 3857245.6653132];    % 20170102_2
% TruePos = [-3026675.29248729 4067187.75522527 3857247.67191648];    % 20170102_1
% TruePos = [-3026675.29412914 4067187.75617734 3857247.6726951];         %rover1
TruePos = [-3026674.48628473 4067187.35187181 3857245.66395674];         %rover2 이동전
% TruePos = [-3026674.57110272 4067187.39555278 3857245.74378219];        % 17년 2월 1일 좌표
% TruePos = [-3026674.33478663 4067187.38883147 3857245.94540826];            % 17년 2월 08일 좌표
% TruePos = [-3026674.33153911 4067187.4137589 3857245.91838126];            % 17년 2월 9일 좌표
% TruePos = [-3026674.33122938 4067187.4134336 3857245.91846153];            % 17년 2월 9일 좌표
% TruePos = [-3026674.33451884 4067187.389255 3857245.94483531];            % 17년 2월 12일 좌표
% [-3026674.57110272,4067187.39555278,3857245.74378219]
% TruePos = [-3026674.28847023 4067187.4147107 3857245.74135094];         %rover2 이동후

% gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2);

%% NMEA text에서 GPGGA, GPGSA, GNGGA 추출
GPGGA = getGPGGA(filename);
GPGSA = getGPGSA(filename);
GN = getGNGGA(filename);

GNGGA_issue = [];
GPGGA_issue = [];
%

%% 결과 분석 시작
gp = 1; gn = 1; count = 0; fixcount = 1;
usedGPS = zeros(length(GPGGA),13);
%% GGA에서 Fix 일때 x,y,z 값 추출해서 TruePos 생성
for i = 1:length(GPGGA)
    
    GGA = GPGGA{i,1};
    %% GPGSA가 존재하면 사용된 GPS 위성 수 저장
    if ~isempty(GPGSA)
        GSA = strsplit(GPGSA{i,1},',');
        
        if str2num(cell2mat(GSA(3))) == 3
            num_gps(i,1) = length(GSA) - 7;
            for s = 4:length(GSA) - 4;
                usedGPS(i,s-2) = str2num(cell2mat(GSA(s)));
            end
        else
            num_gps(i,1) = 0;
        end
    else
        num_gps(i,1) = 0;
    end
    
    [hh,mm,ss,x,y,z,la,lo,qi,nSats,ht] = NEWreadGGA3(GGA) ;
    [gw, gs] = date2gwgs(yy, mo, dd, hh, mm, ss);
    GNGGA(i,1:12) = [round(gs),hh,mm,ss,x,y,z,la,lo,qi,nSats,ht];
    
    if GGA(1:6) == '$GPGGA'
        GNGGA(i,13) = 1;
    elseif GGA(1:6) == '$GNGGA'
        GNGGA(i,13) = 2;
    end
    
    if i> 1 & qi == 4
%     if i> 1 & qi >= 4
        if count > 0
            GNGGA(i-1, 14) = count;
            count = 0;
        end
    elseif qi == 5
        count = count + 1;
    end
    usedGPS(i,1) = round(gs);
end

%% Fixed 결과 sorting
Fixed = GNGGA(find(GNGGA(:,10) >= 4),:);
% Floated = GNGGA(find(GNGGA(:,10) == 5),:);
%% 기준 좌표 계산(Fix Quality 평균값 이용)
% TruePos = [mean(GNGGA(find(GNGGA(:,10)==4),5)), mean(GNGGA(find(GNGGA(:,10)==4),6)), mean(GNGGA(find(GNGGA(:,10)==4),7))];
% TruePos = [mean(Fixed(:,5)), mean(Fixed(:,6)), mean(Fixed(:,7))];
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2);
i=0;cnt=1;
if str2num(filename(5))==2 & str2num(filename(12:14))==112
        Fixed = Fixed([1:6361,6363:length(Fixed)],:);
        disp('Jan, 12, rover2');
    for i=1:length(Fixed)
        if i >= 6362
            TruePos = [-3026674.28847023 4067187.4147107 3857245.74135094];
        else
            TruePos = [-3026674.48628473 4067187.35187181 3857245.66395674];
        end
        gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2);
        GGA = Fixed(i,:);
        dXYZ = [GGA(5), GGA(6), GGA(7)] - TruePos;
        %% dXYZ를 dNEV로 변환
        dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
        %% 각 성분별 RMS 계산
        dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
        dNE = sqrt(dN^2 + dE^2);        %rmsH = myRMS(dNE);
        d2D(i,1) = dNE; d2D(i,2) = Fixed(i,1); d2D(i,3) = GGA(10);
        %rmsV = myRMS(dV);
        d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
        d3D(i,1) = d3; d3D(i,2) = Fixed(i,1); d3D(i,3) = GGA(10);
        result(i,:) = [dN, dE, dV, dNE, d3, GGA(10), num_gps(i,1), GGA(11), Fixed(i,1), i, GGA(13)];
        %     if GNGGA(i,1) == 1 & qi > 4
        %         GPGGA_issue(gp,1) = i;
        %         GPGGA_issue(gp,2:12) = result(i,:);
        %         gp= gp + 1;
        %     elseif GNGGA(i,1) == 2
        %         GNGGA_issue(gn,1) = i;
        %         GNGGA_issue(gn,2:12) = result(i,:);
        %         gn = gn + 1;
        %     end
        cnt=cnt+1;
    end
else
    disp(filename);
    for i=1:length(Fixed)
        GGA = Fixed(i,:);
        dXYZ = [GGA(5), GGA(6), GGA(7)] - TruePos;
        %% dXYZ를 dNEV로 변환
        dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
        %% 각 성분별 RMS 계산
        dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
        dNE = sqrt(dN^2 + dE^2);        %rmsH = myRMS(dNE);
        d2D(i,1) = dNE; d2D(i,2) = Fixed(i,1); d2D(i,3) = GGA(10);
        %rmsV = myRMS(dV);
        d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
        d3D(i,1) = d3; d3D(i,2) = Fixed(i,1); d3D(i,3) = GGA(10);
        result(i,:) = [dN, dE, dV, dNE, d3, GGA(10), num_gps(i,1), GGA(11), Fixed(i,1), i, GGA(13)];
        %     if GNGGA(i,1) == 1 & qi > 4
        %         GPGGA_issue(gp,1) = i;
        %         GPGGA_issue(gp,2:12) = result(i,:);
        %         gp= gp + 1;
        %     elseif GNGGA(i,1) == 2
        %         GNGGA_issue(gn,1) = i;
        %         GNGGA_issue(gn,2:12) = result(i,:);
        %         gn = gn + 1;
        %     end
    end
end

%% Fixed 결과 sorting
result_fixed = result(find(result(:,6) >= 4), :);



d2D_all = rms(d2D(:,1));
d3D_all = rms(d3D(:,1));
d2D_Fix = rms(d2D(find(d2D(:,3) == 4), 1));
d3D_Fix = rms(d3D(find(d3D(:,3) == 4), 1));
d2D_Float = rms(d2D(find(d2D(:,3) == 5), 1));
d3D_Float = rms(d3D(find(d3D(:,3) == 5), 1));
d2D_NF = rms(d2D(find(d2D(:,3) == 1), 1));
d3D_NF = rms(d3D(find(d3D(:,3) == 1), 1));
d2D_all_std = std(d2D(:,1));
d3D_all_std = std(d3D(:,1));
d2D_Fix_std = std(d2D(find(d2D(:,3) == 4), 1));
d3D_Fix_std = std(d3D(find(d3D(:,3) == 4), 1));
d2D_Float_std = std(d2D(find(d2D(:,3) == 5), 1));
d3D_Float_std = std(d3D(find(d3D(:,3) == 5), 1));
d2D_NF_std = std(d2D(find(d2D(:,3) == 1), 1));
d3D_NF_std = std(d3D(find(d3D(:,3) == 1), 1));
% error_result = [decimal((d2D_all)*100)/100, decimal((d3D_all)*100)/100;...
%     decimal((d2D_GPGGA)*100)/100, decimal((d3D_GPGGA)*100)/100;...
%     decimal((d2D_GNGGA)*100)/100, decimal((d3D_GNGGA)*100)/100];
disp(['d2D all = ', num2str(decimal((d2D_all)*100)/100),'m ', '(',num2str(decimal((d2D_all_std)*100)/100),')'])
disp(['d3D all = ', num2str(decimal((d3D_all)*100)/100),'m ', '(',num2str(decimal((d3D_all_std)*100)/100),')'])
disp(['d2D Fix = ', num2str(decimal((d2D_Fix)*100)/100),'m ', '(',num2str(decimal((d2D_Fix_std)*100)/100),')'])
disp(['d3D Fix = ', num2str(decimal((d3D_Fix)*100)/100),'m ', '(',num2str(decimal((d3D_Fix_std)*100)/100),')'])
disp(['d2D Float = ', num2str(decimal((d2D_Float)*100)/100),'m ', '(',num2str(decimal((d2D_Float_std)*100)/100),')'])
disp(['d3D Float = ', num2str(decimal((d3D_Float)*100)/100),'m ', '(',num2str(decimal((d3D_Float_std)*100)/100),')'])
disp(['d2D NF = ', num2str(decimal((d2D_NF)*100)/100),'m ', '(',num2str(decimal((d2D_NF_std)*100)/100),')'])
disp(['d3D NF = ', num2str(decimal((d3D_NF)*100)/100),'m ', '(',num2str(decimal((d3D_NF_std)*100)/100),')'])
disp(['Whole Epoch =', num2str(length(GNGGA(:,1)))]);
disp(['Fixed Epoch =', num2str(length(Fixed(:,1)))]);
disp(['numofFloat =', num2str(length(find(GNGGA(:,14) >0)))])
disp(['meanFloat =', num2str(mean(GNGGA(find(GNGGA(:,14)>5),14)))])
toc;
% RTKlite_postError(GNGGA, result);
% RTKlite_FixedError(Fixed, result_fixed);
fixed_tmp = (['rv',filename(5),filename(11:14),'=Fixed;']);
fixed_result_tmp = (['rv',filename(5),filename(11:14),'result','=result_fixed;']);
eval(fixed_tmp);
eval(fixed_result_tmp);
