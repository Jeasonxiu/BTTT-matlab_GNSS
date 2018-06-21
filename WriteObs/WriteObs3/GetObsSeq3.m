function [ObsSeq] = GetObsSeq3(obs_file,SatType)

% function [ObsSeq, ST, T1] = GetObsSeq3(obs_file)
%
% ������ �ӽ� ����
% ObsSeq: ��(ObsType,10) ��(SatType-G,B,R;3), ������ ������ �� ���ϴ� �����͵��� ��ġ�� �ִ� ���� ��ġ����
% ST(char): ���� ObsType�� G,R,C �κ�
% T1(char): ���� ObsType�� G,R,C �κ��� ���������� �κи�

% Copyright: Mi-So, Kim , January 20th, 2015
% --- Modifications, Hyunu, Tae, January 15, 2016
% �������� : GAL/QZS/SBS �߰�, ������ ������ ������� ������ ����
%% �Լ������ �� test input ����
% clc; clear all;
% obs_file = 'ftna2440.15o';
% [SatType] = GetSatType(obs_file);

fid_obs = fopen(obs_file,'r');
ready = 0; n = 0;
GObs = char('C1C','C1P','C2P','C2C','L1C','L2C','D1C','D2C','S1C','S2C','G');
RObs = char('C1C','C1P','C2P','C2C','L1C','L2C','D1C','D2C','S1C','S2C','R');
CObs = char('C1I','C1Q','C7Q','C7I','L1I','L2Q','D1I','D7I','S1I','S2Q','C');
JObs = char('C1C','C1P','C2P','C2C','L1C','L2C','D1C','D2C','S1C','S2C','J');
EObs = char('C1X','C1P','C2P','C5X','L1X','L5X','D1X','D5X','S1X','S5X','E');
SObs = char('C1C','C1P','C2P','C2C','L1C','L2C','D1C','D2C','S1C','S2C','S');
%%
if ~isempty(SatType)
    if SatType(1) == GObs(11); Obs = GObs;
    elseif SatType(1) == RObs(11); Obs = RObs;
    elseif SatType(1) == CObs(11); Obs = CObs;
    elseif SatType(1) == JObs(11); Obs = JObs;
    elseif SatType(1) == EObs(11); Obs = EObs;
    elseif SatType(1) == SObs(11); Obs = SObs;
    end
    if length(SatType) > 1
        if SatType(2) == GObs(11); Obs = [Obs GObs];
        elseif SatType(2) == RObs(11); Obs = [Obs RObs];
        elseif SatType(2) == CObs(11); Obs = [Obs CObs];
        elseif SatType(2) == JObs(11); Obs = [Obs JObs];
        elseif SatType(2) == EObs(11); Obs = [Obs EObs];
        elseif SatType(2) == SObs(11); Obs = [Obs SObs];
        end
    end
    if length(SatType) > 2
        if SatType(3) == GObs(11); Obs = [Obs GObs];
        elseif SatType(3) == RObs(11); Obs = [Obs RObs];
        elseif SatType(3) == CObs(11); Obs = [Obs CObs];
        elseif SatType(3) == JObs(11); Obs = [Obs JObs];
        elseif SatType(3) == EObs(11); Obs = [Obs EObs];
        elseif SatType(3) == SObs(11); Obs = [Obs SObs];
        end
    end
    if length(SatType) > 3
        if SatType(4) == GObs(11); Obs = [Obs GObs];
        elseif SatType(4) == RObs(11); Obs = [Obs RObs];
        elseif SatType(4) == CObs(11); Obs = [Obs CObs];
        elseif SatType(4) == JObs(11); Obs = [Obs JObs];
        elseif SatType(4) == EObs(11); Obs = [Obs EObs];
        elseif SatType(4) == SObs(11); Obs = [Obs SObs];
        end
    end
    if length(SatType) > 4
        if SatType(5) == GObs(11); Obs = [Obs GObs];
        elseif SatType(5) == RObs(11); Obs = [Obs RObs];
        elseif SatType(5) == CObs(11); Obs = [Obs CObs];
        elseif SatType(5) == JObs(11); Obs = [Obs JObs];
        elseif SatType(5) == EObs(11); Obs = [Obs EObs];
        elseif SatType(5) == SObs(11); Obs = [Obs SObs];
        end
    end
    if length(SatType) > 5
        if SatType(6) == GObs(11); Obs = [Obs GObs];
        elseif SatType(6) == RObs(11); Obs = [Obs RObs];
        elseif SatType(6) == CObs(11); Obs = [Obs CObs];
        elseif SatType(6) == JObs(11); Obs = [Obs JObs];
        elseif SatType(6) == EObs(11); Obs = [Obs EObs];
        elseif SatType(6) == SObs(11); Obs = [Obs SObs];
        end
    end
end
%% ������ ObsType ����
while ~ ready
    s = fgetl(fid_obs);
    if length(s) > 72
        if s(61:73) == 'SYS / # / OBS'
            n = n + 1;
            TypeNum = str2num(s(5:6));
            if TypeNum > 13 % ��ϵ� �������� 13�� �̻��� ��
                T(n,1:58) = s(1:58);
                qq = fix(TypeNum/13);       
                for q = 1 : qq
                s = fgetl(fid_obs);
                T(n,(52*q+8):(52*q+58)) = s(8:58);
                end
            else
                T(n,1:58) = s(1:58);
                %                 T(n,59:109) = ' ';
            end
        end
        if s(61:73) == 'END OF HEADER'
            ready = 1;
            %             fprintf('Error::: Not Found');
        end
    end
end
%% GPS, BDS, GLO, QZSS ObsData ����
for i1 = 1:length(SatType)
    Sa = find(T(:,1) == SatType(i1,1));
    ST(i1,:) = T(Sa,:);
end
Tnum = str2num(ST(:,5:6));
for j1 = 1 : length(Tnum) % ��ġ
    r1 = 1 + 4*(j1-1);
    r2 = r1 + 2;
    for j2 = 1 : Tnum(j1)
        indx1 = 8 + 4*(j2-1);
        indx2 = indx1 + 2;
        T1(j2,r1:r2) = ST(j1,indx1:indx2);
    end
end
%% GPS, BDS, GLO, QZSS ObsData ��ġ ����
ObsNum = length(GObs)-1;
for k1 = 1 : length(Tnum) % 4
%     k1
    for k2 = 1 : ObsNum % 10
%         k2
        g1 = find(T1(:,4*k1-3) == Obs(k2,3*k1-2));
        g2 = find(T1(:,4*k1-2) == Obs(k2,3*k1-1));
        g3 = find(T1(:,4*k1-1) == Obs(k2,3*k1));
        G1 = intersect(g1,g2);
        IN2 = intersect(G1,g3);
        if isempty(IN2)
            continue;
        else
            ObsSeq(k2,k1) = IN2;
        end
    end
end
fclose(fid_obs);