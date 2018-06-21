function [Sat_ar1] = GetSatPosGLO_my(EphGlo,St,deltat)
%
%    <input>
%         eph_glo : GLONASS BRDC array
%         St : 위성위치 계산 시작 시간
%         daltat : 적분간격
%
%    <output>
%        Sat_Array : daltat 간격으로 저장된 Array [Prn-시간-위치-속도-섭동가속도]

%  x = initial value  ex) [1prn 2gs 3x 4y 5z 6xdot 7ydot 8zdot 9xLS 10yLS 11zLS ...]

% - Coded by Mi-So Kim, April 6th 2015

%% 함수만들기 전 테스트
% GloBrdc = 'brdc2760.14g';
% EphGlo = ReadEPH_GLO(GloBrdc);
% St = 438301; deltat = 1;

%%
% Sprn = [1:24]';
Sprn = unique(EphGlo(:,1));
% Sprn = prns;
Sat_ar = zeros(length(Sprn)*1800,11);
intv = 900; kOut = 0;
for i1 = 1:length(Sprn)
    prn = Sprn(i1);
    icol = PickEPH_GLO2(EphGlo,prn,St);
    x = EphGlo(icol,1:11); gs = EphGlo(icol,2);
    x(3:5) = PZ2WGS_JSG(x(3:5));
    x(3:5) = PZ9011toWGS(x(3:5));
    
    dif = St - gs;
    if dif < 0 % 백워드 먼저 → 포워드, 900+900, 최대 1800sec
        it1 = fix(dif/deltat); % integration times: 백워드 적분 횟수1
        it2 = fix(intv/deltat); % : 포워드 적분 횟수
        deltat = -deltat;    sPT = gs;
        % hat = rem(dif,deltat); % 어떻게 적용할지 좀 더 고민
        for k1 = 1:abs(it1) % 백워드
            kOut = kOut + 1;
            sPT = sPT + deltat;
            [xyz,vxyz,xyzLS] = RK4_my(x, deltat);
            Sat_ar(kOut,1) = prn; Sat_ar(kOut,2) = sPT;
            Sat_ar(kOut,3:5) = xyz; Sat_ar(kOut,6:8) = vxyz; Sat_ar(kOut,9:11) = xyzLS;
            x = Sat_ar(kOut,:);
        end
        deltat = -deltat;  sPT = gs;
        x = EphGlo(icol,1:11);
        kOut = kOut + 1;
        Sat_ar(kOut,:) = x;
        for k2 = 1:it2 % 포워드
            kOut = kOut + 1;
            sPT = sPT + deltat;
            [xyz,vxyz,xyzLS] = RK4_my(x, deltat);
            Sat_ar(kOut,1) = prn; Sat_ar(kOut,2) = sPT;
            Sat_ar(kOut,3:5) = xyz; Sat_ar(kOut,6:8) = vxyz; Sat_ar(kOut,9:11) = xyzLS;
            x = Sat_ar(kOut,:);
        end
    else % 포워드 만, 최대 900sec
        it3 = fix(dif/deltat); sPT = gs;
        for k3 = 1:900
            %         for k3 = 1:it3
            kOut = kOut + 1;
            sPT = sPT + deltat;
            [xyz,vxyz,xyzLS] = RK4_my(x, deltat);
            Sat_ar(kOut,1) = prn; Sat_ar(kOut,2) = sPT;
            Sat_ar(kOut,3:5) = xyz; Sat_ar(kOut,6:8) = vxyz; Sat_ar(kOut,9:11) = xyzLS;
            x = Sat_ar(kOut,:);
        end
    end
end

nZeros = find(Sat_ar(:,1) ~= 0);
Sat_ar0 = Sat_ar(nZeros,:);
Sat_ar1 = sortrows(Sat_ar0,2);

% 간단한 확인 작업
% Gprn = 1;
% GLget = find(Sat_ar(:,1) == Gprn); GSst = Sat_ar(GLget,:);
% AA = sortrows(GSst,2);
%
% figure(prn)
% subplot(3,1,1)
% plot(AA(:,2),AA(:,3));
% subplot(3,1,2)
% plot(AA(:,2),AA(:,4));
% subplot(3,1,3)
% plot(AA(:,2),AA(:,5));