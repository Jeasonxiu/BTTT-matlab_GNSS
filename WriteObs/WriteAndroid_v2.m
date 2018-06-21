% function [qm] = WriteAndroid2(filename)
%
%   function [] = WriteAndroid(filename)
%   Read the given logged Android file and make a 'QMfile'
%
%   ex) WriteAndroid('pseudoranges_log_2017_02_14_18_58_59.txt');
%
%   Coded by Joonseong Gim, Feb 16, 2017
%
%   �߰� ���� : Raw Measurement �� State���� ���� ŭ���� �ұ��ϰ� PseudoRange ���� �̻��� ������
%   �����Ǵ°��� �߰���.
%   �׷��� ù Epoch���� State �� �� �ִ밪�� ã�� ��, States ����� 3���� ReceivedSvTimeNanos ����
%   �߰��ϰ�, ReceivedSvTimeNanos ���� 1,000,000,000,000 �̻�ɶ�, State�� �� �ִ밪�� ã����
%   ������.

clear all
close all
% filename = 'pseudoranges_log_2017_02_14_18_58_59.txt';
% filename = 'pseudoranges_log_2017_03_02_15_33_59.txt';      % �������� �޾��ֽ� �ڷ�
% filename = 'pseudoranges_log_2017_08_01_15_39_28.txt';      % �뼺 B ���� ����
% filename = 'pseudoranges_log_2017_08_01_15_50_15.txt';      % �뼺 B ���� �����
% filename = 'pseudoranges_log_2017_08_01_16_02_35.txt';      % �뼺 B ���� ���� ����
filename = 'pseudoranges_log_2017_08_01_16_14_00.txt';      % �뼺 B ���� ���� ����

%% ���� ���
CCC = 299792458;                                    % Speed of Light
WeekSecond = 604800;                                % Week Second
%% Type ���� ���
C1 = 20;                                                                 % Code
C1PrSigmaM = 21;                                                         % PrSigmaM
C1prrSigmaMps = 22;                                                      % PseudorangeRateUncertaintyMetersPerSecond
L1 = 11;                                                                 % reserved
S1 = 41;                                                                 % SNR
U1 = 51;                                                                 % Uncertainty

%% text ���� ���� �� QMfile ����
fid = fopen(filename,'r');
fid_out = fopen('QMfile','w');

%% header ���� ����
ready = 0;
while ready == 0
    s = fgets(fid);
    if length(s) > 100
        if ~isempty(strfind(s,'ElapsedRealtimeMillis'))
            ready = 1;
        end
    end
end
%% QMfile ���� ���� ����
cnt = 0; cnt2=1; Epoch = 0; Line = 1;
States=[];
while 1
    line = fgets(fid);
    gga = {};
    %    while 1
    %% line�� ���� �� ���õ��������� Ȯ���ؼ� cell�� ����� ����
    if length(line) >= 72
        if line(1:3) == 'Raw'
            cnt = cnt + 1;
            gga = strsplit(line,',');       % GGA cell array
            GGA = textscan(line(1:end),'%s%f%d64%f%f%d64%f%f%f%f%f%f%f%f%d64%d64%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
            %% ������� cell�� ������ �� �׸����� �и��ϴ� ����
            TimeNanos = GGA{3};
            TimeOffsetNanos =  GGA{13};
            FullBiasNanos = GGA{6};
            BiasNanos = GGA{7};
            ReceivedSvTimeNanos = GGA{15};
            ReceivedSvTimeUncertaintyNanos = GGA{16};
            State = GGA{14};
            %% �׹��ý��� ������ ���� Indicator
            %% 1=GPS, 2=SBAS, 3=GLONASS, 4=QZSS, 5=BeiDou, 6=Galileo, 0=Unknown
            ConstellationNum = GGA{29};                                             % Constellation Num
            %% �׹��ý��ۿ� ���� QMfile type �� ���� ����
            if ConstellationNum == 1
                TYPE = 100;                                                                     % GPS
            elseif ConstellationNum == 5
                TYPE = 200;                                                                     % BDS
            elseif ConstellationNum == 3
                TYPE = 300;                                                                     % GLO
            else
                TYPE = 600;                                                                     % ETC
            end
        elseif line(1:3) == 'Fix'
            Epoch = Epoch + 1;
        end
    end
    if  Epoch == 1 & ~isempty(gga)
        States(Line,1:3) = [State, ConstellationNum,ReceivedSvTimeNanos];
        Line = Line + 1;
    end
    %     if line(1:3) == 'fix';
    %         break;end
    %    end
    %% Raw Measurement ���� �Ǵ� ����
    % �Ķ������, State ���� �ִ��� ������ ��ȿ�� ������
    % �߰� : 2017�� 08�� 01�� ���������Ϳ��� States ���� ���� ŭ���� �ұ��ϰ� ReceivedSvTimeNanos��
    % �̻��� ���� �߰�,
    % ���� �ý����� State �ִ밪�� �ٸ��Ƿ� ���� �׹� �ý����� ��� ���� �Ǵ��ؾ���
    if ~isempty(States)
        States = States(find(States(:,3) > 100000000000),:);
        Systems = unique(States(:,2));
        Sys = zeros(length(Systems),2);
        Sys(:,1) =Systems;
        for j = 1:length(Systems)
            Sys(j,2) = max(States(find(States(:,2) == Systems(j)),1));
        end
    end
    
    %% tRxNanos�� ����ϱ� ���� ù epoch�� FullBaisNanos ����
    if cnt == 1
        FullBiasNanosFirstEpoch = int64(FullBiasNanos);
    end
    %% �ι�° Epoch���� PseudoRange ���
    if Epoch > 1
        if ~isempty(gga) &  ReceivedSvTimeUncertaintyNanos < 500 & State == Sys(find(Sys(:,1) == ConstellationNum),2)
            %% ���õ����͸� numeric type�� ���ߴ� ����
            TimeNanos = GGA{3};                                                         % TimaNanos
            TimeOffsetNanos =  GGA{13};                                                 % TimeOffsetNanos
            FullBiasNanos = GGA{6};                                                     % FullBiasNanos
            BiasNanos = GGA{7};                                                         % BiasNanos
            ReceivedSvTimeNanos = GGA{15};                                              % ReceivedSvTimeNanos
            ReceivedSvTimeUncertaintyNanos = GGA{16};                                   % ReceivedSvTimeUncertaintyNanos
            PrSigmaM = double(ReceivedSvTimeUncertaintyNanos)*1e-9*CCC;                 % PrSigmaM
            PseudorangeRateUncertaintyMetersPerSecond = double(GGA{19});                % PseudorangeRateUncertaintyMetersPerSecond
            CNo = GGA{17};                                                              % C/No(SNR)
            PRN = GGA{12};                                                              % PRN(SvID)
            AccumulatedDeltaRangeMeters = GGA{21};
            %% GLONASS���� FCN�� ���
            if ConstellationNum == 3 & PRN > 24
                PRN = PRN - 92;
            end
            WeekNumber = floor(-double(FullBiasNanos)*1e-9/WeekSecond);                 % WeekNumber
            WeekNumberNanos = int64(WeekNumber)*int64(WeekSecond*1e9);                  % WeekNumberNanos
            allRxMillisecond = (((TimeNanos-FullBiasNanos)*1e-6));                      % allRxMillisecond
            FctSeconds = double(allRxMillisecond)*1e-3;                                 % Full cycle time tag of M batched measurements
            
            %% Pseudorange ���� ����
            %% compute tRxNanos using gnssRaw.FullBiasNanos(1), so that tRxNanos includes rx clock drift since the first epoch:
            tRxNanos = (TimeNanos+TimeOffsetNanos)-(FullBiasNanosFirstEpoch+BiasNanos)-WeekNumberNanos;     % tRxNanos
            tRxSeconds = double(tRxNanos)*1e-9;                                                             % tRxNanos : nano sec -> sec
            if ConstellationNum == 3
                %% GLONASS �� Time of day�� ����ϸ�, Moscov �� UTC ���� 3�ð� ����, ���� 18�� ���
                Pseudorangeinsecond = tRxSeconds-(double(ReceivedSvTimeNanos)*1e-9+86400*fix(tRxSeconds/86400)...
                    - 3600*3 +18);                              % Pseudorange in sacond(ms)
            elseif ConstellationNum ~= 3
                Pseudorangeinsecond = tRxSeconds-double(ReceivedSvTimeNanos)*1e-9;                              % Pseudorange in sacond(ms)
            end
            %% Week Rollover ���� �ذ��� ���� ����
            iRollover = Pseudorangeinsecond > WeekSecond/2;
            if any(iRollover)
                disp('WARNING: week rollover detected in time tags. Adjusting')
                prS = Pseudorangeinsecond(iRollover);
                delS = round(prS/WeekSecond)*WeekSecond;
                prS = prS - delS;
                maxBiasSecond = 10;
                if any(prS>maxBiasSecond)
                    disp('Failed to correct week rollover')
                else
                    Pseudorangeinsecond(iRollover) = prS;
                    tRxNanos(iRollover) = tRxNanos(iRollover) - delS;
                    disp('Corrected week rollover')
                end
            end
            %% Pseudorange �� CarreirPhaseCounter ���� ����
            Pseudorange = Pseudorangeinsecond*CCC;                                              % Pseudorange(m)
            CarreirPhaseCount = 0;                                                              % reserved
            
            %% QMfile ����
%             fprintf(fid_out,'%8.3f %4d %5d %16.3f \n',tRxSeconds, PRN + TYPE, C1 + TYPE, Pseudorange);
%             fprintf(fid_out,'%8.3f %4d %5d %16.3f \n',tRxSeconds, PRN + TYPE, L1 + TYPE, AccumulatedDeltaRangeMeters);
%             fprintf(fid_out,'%8.3f %4d %5d %16.3f \n',tRxSeconds, PRN + TYPE, S1 + TYPE, CNo);
            fprintf(fid_out,'%8.3f %4d %5d %16.3f \n',tRxSeconds, PRN, C1 + TYPE, Pseudorange);
            fprintf(fid_out,'%8.3f %4d %5d %16.3f \n',tRxSeconds, PRN, L1 + TYPE, AccumulatedDeltaRangeMeters);
            fprintf(fid_out,'%8.3f %4d %5d %16.3f \n',tRxSeconds, PRN, S1 + TYPE, CNo);
            fprintf(fid_out,'%8.3f %4d %5d %16.8f \n',tRxSeconds, PRN, C1PrSigmaM + TYPE, PrSigmaM);
            fprintf(fid_out,'%8.3f %4d %5d %16.10f \n',tRxSeconds, PRN, C1prrSigmaMps + TYPE, PseudorangeRateUncertaintyMetersPerSecond);
            qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, C1 + TYPE, Pseudorange];
            cnt2=cnt2+1;
            qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, S1 + TYPE, CNo];
            cnt2=cnt2+1;
            %         qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, U1 + TYPE, double(PrSigmaM)];
            qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, L1 + TYPE, AccumulatedDeltaRangeMeters];
            cnt2=cnt2+1;
            qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, C1PrSigmaM + TYPE, PrSigmaM];
            cnt2=cnt2+1;
            qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, C1prrSigmaMps + TYPE, PseudorangeRateUncertaintyMetersPerSecond];
            cnt2=cnt2+1;
        end
        if line == -1, break;  end
    end
end
fclose('all');