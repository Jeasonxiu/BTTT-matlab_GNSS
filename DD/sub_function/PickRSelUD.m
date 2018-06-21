function [SatsEl, highELprn] = PickRSelUD(gs, Sats, eph, vec_site)
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
% Copyright: Kwan-Dong Park @PPSoln, 3/21/2018 
%--- Modifications ----
% 11/8/14 함수이름에 RS(기준위성)이 포함되도록 변경함
% 11/23/14 위성 고도각을 저장해서 리턴하도록 함. 임계고도각 이하 제거 목적임
% 3/21/21 gs, gs+1 시각의 위성 고도각을 비교하여 상승/하강 인덱스를 생성

%% 지평좌표계 기준의 고도각 방위각 계산을 위한 위경도 변환
gd = xyz2gd(vec_site); AppLat = gd(1); AppLon = gd(2);
%% 반복계산 과정에서 갱신할 최대 고도각의 초기치
highEL = 0.;
%% 위성갯수 결정 및 출력 변수 초기화
NoSats = length(Sats);
SatsEl = NaN(NoSats, 3);
SatsEl(:, 1) = gs;
SatsUD = NaN(NoSats, 3);
SatsUD(:, 1) = gs;
%% 반복계산 과정에서 최대고도각 위성 검출
for kS = 1:NoSats;
    prn = Sats(kS);
    % gs 시각의 위성 고도각 계산
    icol = PickEPH(eph, prn, gs);
    vec_sat = GetSatPosNC(eph, icol, gs);
    vec_rho = vec_sat - vec_site;
    [az, el] = xyz2azel(vec_rho, AppLat, AppLon);
    % gs+1 시각의 위성 고도각 계산
    icol_prd = PickEPH(eph, prn, gs+1);
    vec_sat_prd = GetSatPosNC(eph, icol_prd, gs+1);
    vec_rho_prd = vec_sat_prd - vec_site;
    [az_prd, el_prd] = xyz2azel(vec_rho_prd, AppLat, AppLon);
    
    if el > highEL
        highEL = el;
        highELprn = kS;
    end
    if el-el_prd < 0
        UD = 1; 
    else
        UD = 0;
    end
    
    SatsEl(kS, 2) = prn;
    SatsEl(kS, 3) = el;
    SatsEl(kS, 4) = UD;
        
    
%     if el < 15
%         disp([gs prn el])
%     end
%     fprintf(1, '%8d %3d %6.2f \n', gs, prn, el)
end