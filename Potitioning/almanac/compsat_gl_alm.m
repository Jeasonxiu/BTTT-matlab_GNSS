function [p_sat] = compsat_gl_alm(t_GPS, gl_alm, i_alm)

% function [p_sat] = com_p_sat_alm_gl(t_GLO, N, gl_alm,i_alm)
%
% t_UTC: [년 월 일 시 분 초], d_N, gl_alm: 알마낙 배열, i_alm= 계산하려는 알마낙 행 넘버
%
% Computes the XYZ coordinates of a given GPS satellite
%
% Oct 18, 2008
%
% [Hye-In Kim, Inha University]

PRN = gl_alm(i_alm,01); %PRN Number
year = gl_alm(i_alm,02); %year
month = gl_alm(i_alm,03); %month
day = gl_alm(i_alm,04); %day
t_oa = gl_alm(i_alm,05); % the time of obtaining almanac from the 24 hours, with UTC(seconds)
nl = gl_alm(i_alm,06); % frequency slot(-7-24)
health = gl_alm(i_alm,07); % health(0-1)
T_a = gl_alm(i_alm,08); % Equator time(sec)
gl_utc = gl_alm(i_alm,09); % correction GLONASS-UTC(sec)
gp_gl = gl_alm(i_alm,10); % correction GPS-GLONASS(sec)
Clock_offset = gl_alm(i_alm,11); % Correction to Board Time Scale(sec)
Lam_a = gl_alm(i_alm,12); % longitude of the ascending node(radian)
del_i = gl_alm(i_alm,13); % correction of inclination(radian)
omega_a = gl_alm(i_alm,14); % argument of perigee(radian)
e_a = gl_alm(i_alm,15); % eccentricity
del_t = gl_alm(i_alm,16);  % correction to the draconitic period(sec)
del_tt = gl_alm(i_alm,17); % Rate of Draconic period variation(sec)

% % 최근 윤년으로부터의 알마낙 기준시가 포함된 날짜
% N_A = [year, month, day, 0, 0, 0];
% % 최근 윤년으로부터의 계산하고자 하는 시각이 포함된 날짜
% N_A_doy = ymd2doy(N_A(1), N_A(2), N_A(3), N_A(4), N_A(5), N_A(6));
% N_doy = ymd2doy(Calendar(1), Calendar(2), Calendar(3), Calendar(4), Calendar(5), Calendar(6));
% d_N = N_doy-N_A_doy;
% 
% t_UTC = Calendar(4)*3600 + Calendar(5)*60 + Calendar(6);

N_a = CD2WS([year, month, day, 0, 0, 0]);
N = t_GPS;
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

E = Kepler_solve(M,e_a,5);

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
p_sat = p_sat';
