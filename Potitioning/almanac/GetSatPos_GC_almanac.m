function [SatPos] = GetSatPos_GC_almanac(alm, prn, gs)
%
%function [SatPos] = GetSatPos_GC_almanac(alm, prn, gs)
%
%** 알마낙 데이터를 이용해 위성위치 산출
%
%

% clear all
% close all
% 
% load('almanac_17255.mat');
% 
% alm = [alm_gps;alm_bds];
% prn = 101;
% gs = 440000;

e         = alm(find(alm(:,1) == prn),3);
t_oa      = alm(find(alm(:,1) == prn),4);
i_0       = alm(find(alm(:,1) == prn),5);
Omega_dot = alm(find(alm(:,1) == prn),6);
sqtA      = alm(find(alm(:,1) == prn),7);
Omega_0   = alm(find(alm(:,1) == prn),8);
omega     = alm(find(alm(:,1) == prn),9);
M_0       = alm(find(alm(:,1) == prn),10);
Af0       = alm(find(alm(:,1) == prn),11);
Af1       = alm(find(alm(:,1) == prn),12);
week_no   = alm(find(alm(:,1) == prn),13);      % week_no   = eph(i_eph,17); 현재로선 불필요한 부분

if prn < 200
    
    %% WGS84에서 정의된 상수
    mu = 3.986005e14;
    Omega_dot_e = 7.2921151467e-5;
    %% 이후 코드에서 사용하게 될 t_k 정의: 단순히 gs에서 t_oa를 뺌
    t_k = gs - t_oa;
    %% 위성궤도 연산 시작(참고문헌은 다양함, 여기서 굳이 코드 설명을 붙일 필요가 없다고 봄
    A = sqtA^2;
    n_0 = sqrt(mu / (A^3));
    M_k = M_0 + n_0*t_k;
    E_k = ecce_anom(M_k, e, 3);
    nu_k = atan2((sqrt(1 - e^2)*sin(E_k)), (cos(E_k) - e));
    
    p = A*(1-e^2);
    r = p/(1+e*cos(nu_k));
    x = r*cos(nu_k);
    y = r*sin(nu_k);
    
    Omega_k = Omega_0+((Omega_dot-Omega_dot_e).*t_k)-(Omega_dot_e.*t_oa);
    if prn < 200
        T1 = [cos(Omega_k) sin(Omega_k) 0; -sin(Omega_k) cos(Omega_k) 0; 0 0 1];
        T2 = [1 0 0;0 cos(i_0) sin(i_0); 0 -sin(i_0) cos(i_0)];
        T3 = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];
    elseif prn > 200
        if prn > 205
            i = 0.3*pi + i_0;
        else
            i = 0.00*pi + i_0;
        end
        T1 = [cos(Omega_k) sin(Omega_k) 0; -sin(Omega_k) cos(Omega_k) 0; 0 0 1];
        T2 = [1 0 0;0 cos(i) sin(i); 0 -sin(i) cos(i)];
        T3 = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];
    end
    T = T3 * T2 * T1;

    
    [SatPos] = T'*[x; y; 0];
    SatPos = SatPos';
elseif prn > 200 & prn < 300
    % ICD-BDS 에서 정의된 상수
    % Table 5-14 Almanac algorithms for users 참고
    mu = 3.986004418e14;
    Omega_dot_e = 7.2921150e-5;
    A = sqtA^2;
    n_0 = sqrt(mu / (A^3));
    t_k = gs - t_oa;
    M_k = M_0 + n_0*t_k;
    E_k = ecce_anom(M_k, e, 20);
    nu_k = atan2(sqrt(1 - e^2)*sin(E_k)/(1-e*cos(E_k)), (cos(E_k) - e)/(1-e*cos(E_k)));
    Phi_k = nu_k + omega;
    r = A*(1-e*cos(E_k));
    x = r*cos(Phi_k);
    y = r*sin(Phi_k);
    
    if prn > 205
        Omega_k = Omega_0 + (Omega_dot - Omega_dot_e)*t_k - Omega_dot_e*t_oa;
        i = 0.3*pi + i_0;
        x_k = x*cos(Omega_k) - y*cos(i)*sin(Omega_k);
        y_k = x*sin(Omega_k) + y*cos(i)*cos(Omega_k);
        z_k = y*sin(i);
    else
        Omega_k = Omega_0 + Omega_dot*t_k - Omega_dot_e*t_oa;
        i = 0.00*pi + i_0;
        x_gk = x*cos(Omega_k) - y*cos(i)*sin(Omega_k);
        y_gk = x*sin(Omega_k) + y*cos(i)*cos(Omega_k);
        z_gk = y*sin(i);
        Arg_Rx = -5;
        Arg_Rz = Omega_dot_e*t_k;
        Rx = [1 0 0; 0 cosd(Arg_Rx) sind(Arg_Rx); 0 -sind(Arg_Rx) cosd(Arg_Rx)];
        Rz = [cos(Arg_Rz) sin(Arg_Rz) 0; -sin(Arg_Rz) cos(Arg_Rz) 0; 0 0 1];
        
        xyz_k = Rz*Rx*[x_gk; y_gk; z_gk];
        x_k = xyz_k(1);
        y_k = xyz_k(2);
        z_k = xyz_k(3);
    end

    
    [SatPos] = [x_k;y_k;z_k]';
end