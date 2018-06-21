% function [QMnewBs, QMnewRv, QMnew] = DDSkyplot(QM1, QM2, eph, Base, Rover)
%
% function [Bs_azel, Rv_azel] = DDSkyplot(QM1, QM2, Base, Rover)
%
%   Read the QMfiles and then SkyPlot 'Base', 'Rover'
%   
%   input QM1 : Base QMfile
%   input QM2 : Rover QMfile
%
%   output QMnewBs = Base SV's az, el
%   output QMnewRv = Rover SV's az, el
%   output QMnew = common SV's az, el
%
%   Example : [Bs_azel, Rv_azel] = DDSkyplot(QM1, QM2, Base, Rover)
%
%   coded by Joonseong Gim, Mar 9, 2016
%


%% Base, Rover 공통 시간 추출
FinalTTs = intersect(Base(:, 1), Rover(:, 1));
% FinalTTs = [29570:29590,1];

%% Base xyz -> gd변환
for i = 1: length(FinalTTs)
    AppPosBs = Base(i,2:4);
    geodBs(i,:) = xyz2gd(AppPosBs);
    AppLatBs(i,:) = geodBs(i,1);
    AppLonBs(i,:) = geodBs(i,2);
end
%% Rover xyz -> gd변환
for i = 1: length(FinalTTs)
    AppPosRv = Rover(i,2:4);
    geodRv(i,:) = xyz2gd(AppPosRv);
    AppLatRv(i,:) = geodRv(i,1);
    AppLonRv(i,:) = geodRv(i,2);
end


%% 방위각-고도각 계산하고 txt 파일로 저장하기
cntBs = 0; cntRv = 0; cnt = 0;

%% 위성 az, el 계산파트(Base, Rover 각각 계산)
for i = 1: length(FinalTTs)
%     if mod(FinalTTs(i),10) ~= 0
%         continue;
%     end
    gs = FinalTTs(i);
    subQMBs = QM1(find(QM1(:,1)==FinalTTs(i)),:);               % 해당 시각 Base QM Sorting
    subQMRv = QM2(find(QM2(:,1)==FinalTTs(i)),:);               % 해당 시각 Rover QM Sorting
    existprnBs = intersect(unique(eph(:,18)), subQMBs(:,2));
    arrSVBs = zeros(length(existprnBs),1);
    for kk = 1:length(existprnBs)
        arrSVBs(kk) = find(subQMBs(:,2) == existprnBs(kk,1));
    end
    subQMBs = subQMBs(sort(arrSVBs),:);
    existprnRv = intersect(unique(eph(:,18)), subQMRv(:,2));
    arrSVRv = zeros(length(existprnRv),1);
    for kkk = 1:length(existprnRv)
        arrSVRv(kkk) = find(subQMRv(:,2) == existprnRv(kkk,1));
    end
    subQMRv = subQMRv(sort(arrSVRv),:);
    NoSatsBs = length(subQMBs(:,1));                            % 해당 시각 Sorted Base QM number of SV
    NoSatsRv = length(subQMRv(:,1));                            % 해당 시각 Sorted Rover QM number of SV
    Sats = intersect(subQMBs(:, 2), subQMRv(:, 2));               % 해당 시각 Base, Rover 같은 위성 선택
    NoSats = length(Sats(:,1));                                 % 해당 시각 공통 위성 수
    No(i,:) = NoSats;
    %% Base 에서 보이는 위성 az, el
    for k = 1 : NoSatsBs
        prnBs = subQMBs(k,2);
        icolBs = PickEPH(eph, prnBs, gs);           % Pick up the proper column in the ephemerides array
        SatXYZBs = GetSatPosNC(eph, icolBs, gs);    % Compute the XYZ coordinates of a given Base GPS satellite
        dPosBs = SatXYZBs - AppPosBs;
        [azBs,elBs] = xyz2azel(dPosBs,AppLatBs(i),AppLonBs(i)); % Compute Base's Azimuth and Elevation
        cntBs = cntBs + 1;
        QMnewBs(cntBs,:) = [gs, prnBs, azBs, elBs];
    end
    %% Rover 에서 보이는 위성 az, el
    for j = 1 : NoSatsRv
        prnRv = subQMRv(j,2);
        icolRv = PickEPH(eph, prnRv, gs);           % Pick up the proper column in the ephemerides array
        SatXYZRv = GetSatPosNC(eph, icolRv, gs);    % Compute the XYZ coordinates of a given Rover GPS satellite
        dPosRv = SatXYZRv - AppPosRv;
        [azRv,elRv] = xyz2azel(dPosRv,AppLatRv(i),AppLonRv(i)); % Compute Rover's Azimuth and Elevation
        cntRv = cntRv + 1;
        QMnewRv(cntRv,:) = [gs, prnRv, azRv, elRv];
    end
    %% Base, Rover 에서 공통으로 보이는 위성 az, el   
    for s = 1 : NoSats
        prn = Sats(s);
        icol = PickEPH(eph, prn, gs);           % Pick up the proper column in the ephemerides array
        SatXYZ = GetSatPosNC(eph, icol, gs);    % Compute the XYZ coordinates of a given GPS satellite
        dPos = SatXYZ - AppPosBs;
        [az,el] = xyz2azel(dPos,AppLatBs(i),AppLonBs(i)); % Compute Azimuth and Elevation
        cnt = cnt + 1;
        QMnew(cnt,:) = [gs, prn, az, el, No(i,:)];
    end
end
site = ' ';
cutoff = 15;

AzBs = QMnewBs(:,3);     ElBs = QMnewBs(:,4);                   
AzRv = QMnewRv(:,3);     ElRv = QMnewRv(:,4);
% Az = QMnew(:,3);     El = QMnew(:,4);
yyBs = (ElBs-90).* -(cos(AzBs*pi/180));
xxBs = (ElBs-90).* -(sin(AzBs*pi/180));
yyRv = (ElRv-90).* -(cos(AzRv*pi/180));
xxRv = (ElRv-90).* -(sin(AzRv*pi/180));
% yy = (El-90).* -(cos(Az*pi/180));
% xx = (El-90).* -(sin(Az*pi/180));

% figure(33);
% Skymap(cutoff, site);
% 
% figure(33)
% hold on;
% 나타내고자하는 방위각, 고도각의 위치에 빨간색 점을 찍는다.
% plot(xxBs, yyBs, '+','color', 'r');
% plot(xxRv, yyRv, 'd','color', 'b');


figure(44);
Skymap(cutoff, site);    

existprn = unique(QMnew(:,2));
for aa = 1:length(existprn)
    prn = existprn(aa);
    if prn < 10
        PRN = strcat('G0',num2str(prn));
    else
        PRN = strcat('G',num2str(prn));
    end
    indxprn = find(QMnew(:,2) == prn);
    QM = QMnew(indxprn,:);
    Az = QM(:,3);     El = QM(:,4);
    yy = (El-90).* -(cos(Az*pi/180));
    xx = (El-90).* -(sin(Az*pi/180));
    ylast = yy(length(yy));
    xlast = xx(length(xx));
    figure(44);
    hold on;
    plot(xx, yy, '.', 'Markersize', 30);
    plot(xlast, ylast, 'r.', 'Markersize', 30);
    text(xlast, ylast, PRN,'Fontsize',15)
end
% figure(44);
% Skymap(cutoff, site);
% figure(44);
% hold on;
% plot(xx, yy, '.', 'Markersize', 15);
