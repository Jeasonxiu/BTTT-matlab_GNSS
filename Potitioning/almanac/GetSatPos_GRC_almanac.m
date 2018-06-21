function [SatPos] = GetSatPos_GRC_almanac(alm, prn, gs)
%
%function [SatPos] = GetSatPos_GRC_almanac(alm, prn, gs)
%
%** 알마낙 데이터를 이용해 위성위치 산출
%
%
% 
% clear all
% close all
% 
% load('almanac_17255.mat');
% 
% alm = [alm_gps;alm_bds;alm_glo];
% prn = 101;
% gs = 440000;
if prn < 300
    
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
    
    %% GPS 위성 좌표 계산
    if prn <200
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
        
        %% BDS 위성 좌표 계산
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
else
    
    %% GLO 위성 좌표 계산
    PRN = alm(find(alm(:,1) == prn),01); %PRN Number
    year = alm(find(alm(:,1) == prn),02); %year
    month = alm(find(alm(:,1) == prn),03); %month
    day = alm(find(alm(:,1) == prn),04); %day
    t_oa = alm(find(alm(:,1) == prn),05); % the time of obtaining almanac from the 24 hours, with UTC(seconds)
    nl = alm(find(alm(:,1) == prn),06); % frequency slot(-7-24)
    health = alm(find(alm(:,1) == prn),07); % health(0-1)
    T_a = alm(find(alm(:,1) == prn),08); % Equator time(sec)
    gl_utc = alm(find(alm(:,1) == prn),09); % correction GLONASS-UTC(sec)
    gp_gl = alm(find(alm(:,1) == prn),10); % correction GPS-GLONASS(sec)
    Clock_offset = alm(find(alm(:,1) == prn),11); % Correction to Board Time Scale(sec)
    Lam_a = alm(find(alm(:,1) == prn),12); % longitude of the ascending node(radian)
    del_i = alm(find(alm(:,1) == prn),13); % correction of inclination(radian)
    omega_a = alm(find(alm(:,1) == prn),14); % argument of perigee(radian)
    e_a = alm(find(alm(:,1) == prn),15); % eccentricity
    del_t = alm(find(alm(:,1) == prn),16);  % correction to the draconitic period(sec)
    del_tt = alm(find(alm(:,1) == prn),17); % Rate of Draconic period variation(sec)
        
%     N_a = CD2WS([year, month, day, 0, 0, 0]);
    N_a = 623;
    N = gs;
    diff_N = N-N_a;
    
    t_GLO = diff_N + 3*3600;
    
    % 궤도 경사각 평균(radian)
    i_mean = 63*pi/180;
    % Draconian period T의 평균(sec)
    T_mean = 43200;
    i = i_mean + del_i;
    T = T_mean + del_t;
    
    a_E = 6378136 ;  % Earth Semi-major axis
    omega_dot_E = 7.292115 * 10^(-5);   %Earth rotation rate
    mu = 3.9860044 * 10^14;  %Gravitational constant
    
    % t_k = (N-N_a)*86400 + t_GLO-T_a;
    t_k = t_GLO-T_a;
    n = 2*pi/T;
    a = (mu/(n^2))^(1/3);
    d_Lam = -10*((a_E/a)^(7/2))*(pi/(180*86400))*cos(i);
    d_omega = 5*((a_E/a)^(7/2))*(pi/(180*86400))*((5*(cos(i))^2)-1);
    Lam = Lam_a + (d_Lam - omega_dot_E)*t_k;
    omega = omega_a + d_omega*t_k;
    E_phi = 2*atan(tan(omega/2)*sqrt((1-e_a)/(1+e_a)));
    if omega < pi
        del_T = (E_phi-(e_a*sin(E_phi)))/n + 0;
    else
        del_T = (E_phi-(e_a*sin(E_phi)))/n + T;
    end
    M = n * (t_k - del_T);
    
    E = ecce_anom(M,e_a,5);
    
    pos = a*[cos(E) - e_a; (sqrt(1-(e_a)^2))*sin(E); 0];
    vel = (a/(1-e_a*cos(E))) * [-n*sin(E); n*(sqrt(1-(e_a)^2))*cos(E); 0];
    
    e1 = [cos(omega)*cos(Lam) - sin(omega)*sin(Lam)*cos(i);
        cos(omega)*sin(Lam) + sin(omega)*cos(Lam)*cos(i);
        sin(omega)*sin(i)];
    
    e2 = [-sin(omega)*cos(Lam) - cos(omega)*sin(Lam)*cos(i);
        -sin(omega)*sin(Lam) + cos(omega)*cos(Lam)*cos(i);
        cos(omega)*sin(i)];
    
    p_sat = pos(1)*e1 + pos(2)*e2;
    v_sat = vel(1)*e1 + vel(2)*e2 + omega_dot_E*[p_sat(2); -p_sat(1); 0];
    [SatPos] = p_sat';
    
end
