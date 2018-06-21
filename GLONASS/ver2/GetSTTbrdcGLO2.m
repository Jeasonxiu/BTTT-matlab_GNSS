function [STT] = GetSTTbrdcGLO2(Sat_ar,gs,prn,vec_site)

% function [STT] = GetSTTbrdcGLO2(Sat_ar,gs,prn,ver_site)
% 
% <input> Sat_ar : GetSatPos_new�� 1�ʰ������� ������� Sat array
%         gs : time-epoch given in GPS Week Second
%         prn : Satellite number
%         ver_site : Site Position
% 
% <output> STT : signal transmission time [second]
% 
% Code by  Mi-So Kim, April 7th, 2015,

%% ����� �Ű����� ����
CCC = 299792458;
eps = 1.e-10;
maxIter = 10;
%% �ݺ����
stt_0 = 0.075;
for iter = 1:maxIter
   [epochSat,PoS,Vel] = findSatPos(Sat_ar,gs,prn);
   stt = -stt_0;
   [xyz,vxyz,xyzLS] = RK4_new(epochSat,stt);
   vec_sat = xyz; 
   vec_sat = vec_sat';
   com = norm(vec_sat - vec_site);
   
   STT = com/CCC;
    if abs(STT - stt_0) < eps
        return
    end
    stt_0 = STT;
end