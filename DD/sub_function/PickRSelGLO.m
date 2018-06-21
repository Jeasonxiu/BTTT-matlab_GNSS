function [SatsEl, highELprn] = PickRSelGLO(gs, LeapSec, Sats, ephglo, Sat_ar, TauC, vec_site)
%
% DO: Find the index of RS(Reference Satellite) with the highest elevation angle
%
% <input>   gs: GPS Week Second
%           Sats: Satellite PRN array
%           eph; Ephemeris array
%           vec_site: Site Vector
%
% <output>  highELprn: Index of RS - PRN with the highest elevation angle
%           SatsEl: List of PRNs with elevation angle
%
% Copyright: Kwan-Dong Park @Jipyong Space, 11/7/2014
%--- Modifications ----
% 11/8/14 함수이름에 RS(기준위성)이 포함되도록 변경함
% 11/23/14 위성 고도각을 저장해서 리턴하도록 함. 임계고도각 이하 제거 목적임

%% 지평좌표계 기준의 고도각 방위각 계산을 위한 위경도 변환
gd = xyz2gd(vec_site); AppLat = gd(1); AppLon = gd(2);
%% 반복계산 과정에서 갱신할 최대 고도각의 초기치
highEL = 0.;
%% 위성갯수 결정 및 출력 변수 초기화
NoSats = length(Sats);
SatsEl = NaN(NoSats, 3);
SatsEl(:, 1) = gs;
%% 반복계산 과정에서 최대고도각 위성 검출
for kS = 1:NoSats;
    prn = Sats(kS);
    
    tc = gs - LeapSec;
    icol=PickEPH_GLO2(ephglo, prn, tc);
    
    TauN=ephglo(icol,12); GammaN=ephglo(icol,13); %: tau & gamma 시계오차 보정에 사용
    ch_num=ephglo(icol,16); %: channel number 전리층 보정에 사용
    
    % 신호전달시간 계산
    STT = GetSTTbrdcGLO2(Sat_ar,gs,prn,vec_site');
    % LeapSecond & 신호전달 시간을 보정한 위성 위치 산출
    [vec_sat, SatVel] = SatPosLS_STT(Sat_ar,gs,prn,LeapSec,STT,TauC);
    vec_rho = vec_sat - vec_site;
    
    [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
    
    
    if el > highEL
        highEL = el;
        highELprn = kS;
    end
    
    SatsEl(kS, 2) = prn;
    SatsEl(kS, 3) = el;
    
    %     if el < 15
    %         disp([gs prn el])
    %     end
    %     fprintf(1, '%8d %3d %6.2f \n', gs, prn, el)
end