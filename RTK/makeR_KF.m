function [R] = makeR_KF(SatsList_OS, SatsInfo_Bs, SatsInfo_Rv, idx_RS, elecut)
%   �����л� R ���� �Լ� �⺻������ ��ȯ

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
    % Other Sat ���� ����
    el_OS_Bs = SatsInfo_Bs(prn_OS, 7);
    el_OS_Rv = SatsInfo_Rv(prn_OS, 7);
    % ����� Other Sat ���� ���
    el_OS_mean = (el_OS_Bs+el_OS_Rv)/2;
    % ��� ������ ������ ����(�����л�����)
    Noise_1 = ( 1/sind(el_OS_mean) )^(1);
    Noise_2 = ( 1/sind(el_OS_mean) )^(2);
    Noise_3 = ( 1/sind(el_OS_mean) )^(3);
    Noise_4 = ( 1/sind(el_OS_mean) )^(5);
    Noise_5 = ( 1/sind(el_OS_mean) )^(10);
    
    % ������ ���� ������ R����
    if (el_RS_Bs > elecut) && (el_RS_Rv > elecut) && (el_OS_Bs > elecut) && (el_OS_Rv > elecut) %% - ������ �ٲ㰡�鼭 �غ���
        if (el_RS_Bs > elecut3) && (el_RS_Rv > elecut3) && (el_OS_Bs > elecut3) && (el_OS_Rv > elecut3) %% - ������ �ٲ㰡�鼭 �غ���
            R(i,i) = Noise_1;
        elseif (el_RS_Bs > elecut2) && (el_RS_Rv > elecut2) && (el_OS_Bs > elecut2) && (el_OS_Rv > elecut2) %% - ������ �ٲ㰡�鼭 �غ���
            R(i,i) = Noise_2;       % Best (17.11.21)
        elseif (el_RS_Bs > elecut1) && (el_RS_Rv > elecut1) && (el_OS_Bs > elecut1) && (el_OS_Rv > elecut1) %% - ������ �ٲ㰡�鼭 �غ���
            R(i,i) = Noise_3;       % Best (17.11.21)
        else
            R(i,i) = Noise_4;
        end
    else
        R(i,i) = Noise_5;
    end
end
