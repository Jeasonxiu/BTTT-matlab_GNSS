function [R] = makeR_KF(SatsList_OS, SatsInfo_Bs, SatsInfo_Rv, idx_RS, elecut)
%   상현학생 R 생성 함수 기본식으로 변환

% clear all
% close all
% load('maker_test.mat');

elecut1 = 30;
elecut2 = 45;
elecut3 = 60;


el_RS_Bs = SatsInfo_Bs(idx_RS);
el_RS_Rv = SatsInfo_Rv(idx_RS);

for i=1:length(SatsList_OS)
    prn_OS = SatsList_OS(i);
    % Other Sat 고도각 추출
    el_OS_Bs = SatsInfo_Bs(prn_OS, 7);
    el_OS_Rv = SatsInfo_Rv(prn_OS, 7);
    % 추출된 Other Sat 고도각 평균
    el_OS_mean = (el_OS_Bs+el_OS_Rv)/2;
    % 평균 고도각별 노이즈 설정(상현학생기준)
    Noise_1 = ( 1/sind(el_OS_mean) )^(1);
    Noise_2 = ( 1/sind(el_OS_mean) )^(2);
    Noise_3 = ( 1/sind(el_OS_mean) )^(3);
    Noise_4 = ( 1/sind(el_OS_mean) )^(5);
    Noise_5 = ( 1/sind(el_OS_mean) )^(10);
    
    % 고도각에 따른 위성별 R생성
    if (el_RS_Bs > elecut) && (el_RS_Rv > elecut) && (el_OS_Bs > elecut) && (el_OS_Rv > elecut) %% - 각도를 바꿔가면서 해보고
        if (el_RS_Bs > elecut3) && (el_RS_Rv > elecut3) && (el_OS_Bs > elecut3) && (el_OS_Rv > elecut3) %% - 각도를 바꿔가면서 해보고
            R(i,i) = Noise_1;
        elseif (el_RS_Bs > elecut2) && (el_RS_Rv > elecut2) && (el_OS_Bs > elecut2) && (el_OS_Rv > elecut2) %% - 각도를 바꿔가면서 해보고
            R(i,i) = Noise_2;       % Best (17.11.21)
        elseif (el_RS_Bs > elecut1) && (el_RS_Rv > elecut1) && (el_OS_Bs > elecut1) && (el_OS_Rv > elecut1) %% - 각도를 바꿔가면서 해보고
            R(i,i) = Noise_3;       % Best (17.11.21)
        else
            R(i,i) = Noise_4;
        end
    else
        R(i,i) = Noise_5;
    end
end
