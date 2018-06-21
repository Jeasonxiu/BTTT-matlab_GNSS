function [SatsEl, highELprn] = PickRSelGloSP3(gs, Sats, SP3, vec_site)
%
% DO: Find the index of RS(Reference Satellite) with the highest elevation angle
% 
% <input>   gs: GPS Week Second
%           Sats: Satellite PRN array
%           SP3; Ephemeris array
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
    vec_sat = IntpSP3e1(SP3, prn, gs); vec_sat = vec_sat(2:4);
    vec_rho = vec_sat - vec_site;
    [az, el] = xyz2azel(vec_rho, AppLat, AppLon);
    
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