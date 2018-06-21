% function [qm] = WriteAndroid2(filename)
%
%   function [] = WriteAndroid2(filename)
%   Read the given logged Android file and make a 'QMfile'
%
%   ex) WriteAndroid('pseudoranges_log_2017_02_14_18_58_59.txt');
%
%   Coded by Joonseong Gim, June 11, 2018
%
%   추가 사항 : Raw Measurement 중 State값이 가장 큼에도 불구하고 PseudoRange 값이 이상한 위성이
%   생성되는것을 발견함.
%   그래서 86번째 라인 if bitand(State,2^0) &  bitand(State,2^3)를 추가하여 State 값을
%   판단함(첫Bit와 4번째 Bit가 모두 1이상)
%   위성군 Type 변경
%   GLONASS state 판단 인자 수정 bitand(State, 2^7)
%   GLONASS QM변환 가능

clear all
close all

% filename = 'S8_B_gnss_log_2018_06_06_16_57_12.txt';    % 2018-06-06 S8-1
% filename = 'S8_B_gnss_log_2018_06_06_18_08_01.txt';    % 2018-06-06 S8-2
% filename = 'N8_B_gnss_log_2018_06_06_16_56_47.txt';    % 2018-06-06 Note8-1
% filename = 'N8_B_gnss_log_2018_06_06_18_07_59.txt';    % 2018-06-06 Note8-2
% filename = 'S8_B_gnss_log_2018_06_08_14_57_31.txt';    % 2018-06-08 S8-1
% filename = 'S8_B_gnss_log_2018_06_08_16_02_16.txt';    % 2018-06-08 S8-2
% filename = 'N8_B_gnss_log_2018_06_08_14_57_29.txt';    % 2018-06-08 Note8-1
% filename = 'N8_B_gnss_log_2018_06_08_16_02_00.txt';    % 2018-06-08 Note8-2
filename = 'gnss_log_2018_06_13_15_45_53.txt';    % 2018-06-13 S9
%% 고정 상수
CCC = 299792458;                                    % Speed of Light
WeekSecond = 604800;                                % Week Second
%% Type 구분 상수
C1 = 20;                                                                 % Code
C1PrSigmaM = 21;                                                         % PrSigmaM
C1prrSigmaMps = 22;                                                      % PseudorangeRateUncertaintyMetersPerSecond
L1 = 11;                                                                 % reserved
D1 = 31;                                                                 % Dopller
S1 = 41;                                                                 % SNR
U1 = 51;                                                                 % Uncertainty

%% GLONASS Freq Num
FCN=[1;-4;5;6;1;-4;5;6;-6;-7;0;-1;-2;-7;0;-1;4;-3;3;2;4;-3;3;2];
%% text 파일 열기 및 QMfile 선언
fid = fopen(filename,'r');
fid_out = fopen('QMfile','w');
fid_out2 = fopen('QMstate','w');

%% header 삭제 과정
ready = 0;
while ready == 0
    s = fgets(fid);
    if length(s) > 100
        if ~isempty(strfind(s,'ElapsedRealtimeMillis'))
            disp('Delete Header');
            ready = 1;
        end
    end
end
%% QMfile 생성 과정 시작
cnt = 0; cnt2=1; Epoch = 0; FirstState = 0;
while 1
    line = fgets(fid);
    gga = {};
    %    while 1
    %% line을 읽은 후 원시데이터인지 확인해서 cell로 만드는 과정
    if length(line) >= 72
        if line(1:3) == 'Raw'
            cnt = cnt + 1;
            gga = strsplit(line,',');       % GGA cell array
            if length(gga) == 20
                GGA = textscan(line(1:end),'%s%f%d64%f%f%d64%f%f%f%f%f%f%f%f%d64%d64%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
            elseif length(gga) == 22
                GGA = textscan(line(1:end),'%s%f%d64%f%f%d64%f%f%f%f%f%f%f%f%d64%d64%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
            else
                GGA = textscan(line(1:end),'%s%f%d64%f%f%d64%f%f%f%f%f%f%f%f%d64%d64%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
            end
            %% 만들어진 cell을 가지고 각 항목으로 분리하는 과정
            TimeNanos = GGA{3};
            TimeOffsetNanos =  GGA{13};
            FullBiasNanos = GGA{6};
            BiasNanos = GGA{7};
            ReceivedSvTimeNanos = GGA{15};
            ReceivedSvTimeUncertaintyNanos = GGA{16};
            State = GGA{14};
            %% 항법시스템 구분을 위한 Indicator
            %% 1=GPS, 2=SBAS, 3=GLONASS, 4=QZSS, 5=BeiDou, 6=Galileo, 0=Unknown
            ConstellationNum = GGA{29};                                             % Constellation Num
            %% 항법시스템에 따른 QMfile type 값 단위 설정
            switch ConstellationNum
                case 1
                    TYPE = 100;                                                                 % GPS
                    f = 1575.42e6;
                    WL = CCC/f;
                case 2
                    TYPE = 600;                                                                 % SBAS
                    f = 1575.42e6;
                    WL = CCC/f;
                case 3
                    TYPE = 300;                                                                 % GLO
                    % 아래에서 PRN에 맞는 채널 번호로 파장값 산출
                case 4
                    TYPE = 500;                                                                 % QZSS
                    f = 1575.42e6;
                    WL = CCC/f;
                case 5
                    TYPE = 200;                                                                 % BDS
                    f = 1589.74e6;
                    WL = CCC/f;
                case 6
                    TYPE = 400;                                                                 % Galileo
                    f = 1575.42e6;
                    WL = CCC/f;
                otherwise
                    TYPE = 700;                                                                 % ETC
            end
            %% 두번째 Epoch부터 PseudoRange 계산
            if bitand(State,2^0) &  bitand(State,2^3)
                %% tRxNanos를 계산하기 위한 첫 epoch의 FullBaisNanos 저장
                if cnt2 == 1
                    FullBiasNanosFirstEpoch = int64(FullBiasNanos);
                    FirstState = State;
                end
            end
            if ConstellationNum ~= 3
                if bitand(State,2^0) &  bitand(State,2^3) & cnt >= 1
                    %% 원시데이터를 numeric type에 맞추는 과정
                    TimeNanos = GGA{3};                                                         % TimaNanos
                    TimeOffsetNanos =  GGA{13};                                                 % TimeOffsetNanos
                    FullBiasNanos = GGA{6};                                                     % FullBiasNanos
                    BiasNanos = GGA{7};                                                         % BiasNanos
                    ReceivedSvTimeNanos = GGA{15};                                              % ReceivedSvTimeNanos
                    ReceivedSvTimeUncertaintyNanos = GGA{16};                                   % ReceivedSvTimeUncertaintyNanos
                    PrSigmaM = double(ReceivedSvTimeUncertaintyNanos)*1e-9*CCC;                 % PrSigmaM
                    PseudorangeRateMetersPerSecond = double(GGA{18});                           % PseudorangeRateMetersPerSecond
                    PseudorangeRateUncertaintyMetersPerSecond = double(GGA{19});                % PseudorangeRateUncertaintyMetersPerSecond
                    CNo = GGA{17};                                                              % C/No(SNR)
                    PRN = GGA{12};                                                              % PRN(SvID)
                    AccumulatedDeltaRangeMeters = GGA{21};
                    %% GLONASS Doppler 계산을 위한 FCN 설정 및 WL계산 및 QZSS PRN 번호 -100
                    if TYPE == 500
                        PRN = PRN - 100;
                    end
                    Doppler = -PseudorangeRateMetersPerSecond / WL;                              % Doppler Measurement
                    WeekNumber = floor(-double(FullBiasNanos)*1e-9/WeekSecond);                 % WeekNumber
                    WeekNumberNanos = int64(WeekNumber)*int64(WeekSecond*1e9);                  % WeekNumberNanos
                    allRxMillisecond = (((TimeNanos-FullBiasNanos)*1e-6));                      % allRxMillisecond
                    FctSeconds = double(allRxMillisecond)*1e-3;                                 % Full cycle time tag of M batched measurements
                    
                    %% Pseudorange 생성 과정
                    %% compute tRxNanos using gnssRaw.FullBiasNanos(1), so that tRxNanos includes rx clock drift since the first epoch:
                    tRxNanos = (TimeNanos+TimeOffsetNanos)-(FullBiasNanosFirstEpoch+BiasNanos)-WeekNumberNanos;     % tRxNanos
                    tRxSeconds = double(tRxNanos)*1e-9;                                                             % tRxNanos : nano sec -> sec
                    Pseudorangeinsecond = tRxSeconds-double(ReceivedSvTimeNanos)*1e-9;                              % Pseudorange in sacond(ms)
                    %% Week Rollover 문제 해결을 위한 과정
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
                    %% Pseudorange 및 CarreirPhaseCounter 생성 과정
                    Pseudorange = Pseudorangeinsecond*CCC;                                              % Pseudorange(m)
                    CarreirPhaseCount = 0;                                                              % reserved
                    
                    %% QMfile 생성
                    if tRxNanos >= 0
                        %             fprintf(fid_out,'%8.3f %4d %5d %16.3f \n',tRxSeconds, PRN + TYPE, C1 + TYPE, Pseudorange);
                        %             fprintf(fid_out,'%8.3f %4d %5d %16.3f \n',tRxSeconds, PRN + TYPE, L1 + TYPE, AccumulatedDeltaRangeMeters);
                        %             fprintf(fid_out,'%8.3f %4d %5d %16.3f \n',tRxSeconds, PRN + TYPE, S1 + TYPE, CNo);
                        fprintf(fid_out,'%8.9f %4d %5d %16.7f \n',tRxSeconds, PRN + TYPE, C1 + TYPE, Pseudorange);
                        fprintf(fid_out,'%8.9f %4d %5d %16.11f \n',tRxSeconds, PRN + TYPE, L1 + TYPE, AccumulatedDeltaRangeMeters);
                        fprintf(fid_out,'%8.9f %4d %5d %16.3f \n',tRxSeconds, PRN + TYPE, D1 + TYPE, Doppler);
                        fprintf(fid_out,'%8.9f %4d %5d %16.3f \n',tRxSeconds, PRN + TYPE, S1 + TYPE, CNo);
                        fprintf(fid_out,'%8.9f %4d %5d %16.8f \n',tRxSeconds, PRN + TYPE, C1PrSigmaM + TYPE, PrSigmaM);
                        fprintf(fid_out2,'%8.9f %4d %5d\n',tRxSeconds, PRN, GGA{20});
                        %                     fprintf(fid_out,'%8.3f %4d %5d %16.10f \n',tRxSeconds, PRN, C1prrSigmaMps + TYPE, PseudorangeRateUncertaintyMetersPerSecond);
                        qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, C1 + TYPE, Pseudorange];
                        States(cnt2, 1:3) = [tRxSeconds, PRN + TYPE, GGA{20}];
                        cnt2=cnt2+1;
                        qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, S1 + TYPE, CNo];
                        cnt2=cnt2+1;
                        %         qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, U1 + TYPE, double(PrSigmaM)];
                        qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, L1 + TYPE, AccumulatedDeltaRangeMeters];
                        cnt2=cnt2+1;
                        qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, D1 + TYPE, Doppler];
                        cnt2=cnt2+1;
                        qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, C1PrSigmaM + TYPE, PrSigmaM];
                        cnt2=cnt2+1;
                        %                     qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, C1prrSigmaMps + TYPE, PseudorangeRateUncertaintyMetersPerSecond];
                        %                     cnt2=cnt2+1;
                    end
                    
                end
            else
                if bitand(State,2^7) & cnt >= 1
                    %% 원시데이터를 numeric type에 맞추는 과정
                    TimeNanos = GGA{3};                                                         % TimaNanos
                    TimeOffsetNanos =  GGA{13};                                                 % TimeOffsetNanos
                    FullBiasNanos = GGA{6};                                                     % FullBiasNanos
                    BiasNanos = GGA{7};                                                         % BiasNanos
                    ReceivedSvTimeNanos = GGA{15};                                              % ReceivedSvTimeNanos
                    ReceivedSvTimeUncertaintyNanos = GGA{16};                                   % ReceivedSvTimeUncertaintyNanos
                    PrSigmaM = double(ReceivedSvTimeUncertaintyNanos)*1e-9*CCC;                 % PrSigmaM
                    PseudorangeRateMetersPerSecond = double(GGA{18});                           % PseudorangeRateMetersPerSecond
                    PseudorangeRateUncertaintyMetersPerSecond = double(GGA{19});                % PseudorangeRateUncertaintyMetersPerSecond
                    CNo = GGA{17};                                                              % C/No(SNR)
                    PRN = GGA{12};                                                              % PRN(SvID)
                    AccumulatedDeltaRangeMeters = GGA{21};
                    %% GLONASS에서 FCN인 경우
                    if ConstellationNum == 3 & PRN > 24
                        PRN = PRN - 92;
                    end
                    %% GLONASS Doppler 계산을 위한 FCN 설정 및 WL계산 및 QZSS PRN 번호 -100
                    if TYPE == 300
                        f = (1602 + FCN(PRN) * 0.5625) * 1e6;
                        WL = CCC / f;
                    end
                    Doppler = -PseudorangeRateMetersPerSecond / WL;                              % Doppler Measurement
                    WeekNumber = floor(-double(FullBiasNanos)*1e-9/WeekSecond);                 % WeekNumber
                    WeekNumberNanos = int64(WeekNumber)*int64(WeekSecond*1e9);                  % WeekNumberNanos
                    allRxMillisecond = (((TimeNanos-FullBiasNanos)*1e-6));                      % allRxMillisecond
                    FctSeconds = double(allRxMillisecond)*1e-3;                                 % Full cycle time tag of M batched measurements
                    
                    %% Pseudorange 생성 과정
                    %% compute tRxNanos using gnssRaw.FullBiasNanos(1), so that tRxNanos includes rx clock drift since the first epoch:
                    tRxNanos = (TimeNanos+TimeOffsetNanos)-(FullBiasNanosFirstEpoch+BiasNanos)-WeekNumberNanos;     % tRxNanos
                    tRxSeconds = double(tRxNanos)*1e-9;                                                             % tRxNanos : nano sec -> sec
                    %% GLONASS 는 Time of day를 사용하며, Moscov 는 UTC 기준 3시간 빠름, 윤초 18초 고려
                    Pseudorangeinsecond = tRxSeconds-(double(ReceivedSvTimeNanos)*1e-9+86400*fix(tRxSeconds/86400)...
                        - 3600*3 +18);
                    %% Week Rollover 문제 해결을 위한 과정
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
                    %% Pseudorange 및 CarreirPhaseCounter 생성 과정
                    Pseudorange = Pseudorangeinsecond*CCC;                                              % Pseudorange(m)
                    CarreirPhaseCount = 0;                                                              % reserved
                    
                    %% QMfile 생성
                    if tRxNanos >= 0
                        %             fprintf(fid_out,'%8.3f %4d %5d %16.3f \n',tRxSeconds, PRN + TYPE, C1 + TYPE, Pseudorange);
                        %             fprintf(fid_out,'%8.3f %4d %5d %16.3f \n',tRxSeconds, PRN + TYPE, L1 + TYPE, AccumulatedDeltaRangeMeters);
                        %             fprintf(fid_out,'%8.3f %4d %5d %16.3f \n',tRxSeconds, PRN + TYPE, S1 + TYPE, CNo);
                        fprintf(fid_out,'%8.9f %4d %5d %16.7f \n',tRxSeconds, PRN + TYPE, C1 + TYPE, Pseudorange);
                        fprintf(fid_out,'%8.9f %4d %5d %16.11f \n',tRxSeconds, PRN + TYPE, L1 + TYPE, AccumulatedDeltaRangeMeters);
                        fprintf(fid_out,'%8.9f %4d %5d %16.3f \n',tRxSeconds, PRN + TYPE, D1 + TYPE, Doppler);
                        fprintf(fid_out,'%8.9f %4d %5d %16.3f \n',tRxSeconds, PRN + TYPE, S1 + TYPE, CNo);
                        fprintf(fid_out,'%8.9f %4d %5d %16.8f \n',tRxSeconds, PRN + TYPE, C1PrSigmaM + TYPE, PrSigmaM);
                        fprintf(fid_out2,'%8.9f %4d %5d\n',tRxSeconds, PRN, GGA{20});
                        %                     fprintf(fid_out,'%8.3f %4d %5d %16.10f \n',tRxSeconds, PRN, C1prrSigmaMps + TYPE, PseudorangeRateUncertaintyMetersPerSecond);
                        qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, C1 + TYPE, Pseudorange];
                        States(cnt2, 1:3) = [tRxSeconds, PRN + TYPE, GGA{20}];
                        cnt2=cnt2+1;
                        qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, S1 + TYPE, CNo];
                        cnt2=cnt2+1;
                        %         qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, U1 + TYPE, double(PrSigmaM)];
                        qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, L1 + TYPE, AccumulatedDeltaRangeMeters];
                        cnt2=cnt2+1;
                        qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, D1 + TYPE, Doppler];
                        cnt2=cnt2+1;
                        qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, C1PrSigmaM + TYPE, PrSigmaM];
                        cnt2=cnt2+1;
                        %                     qm(cnt2,1:4) = [tRxSeconds, PRN + TYPE, C1prrSigmaMps + TYPE, PseudorangeRateUncertaintyMetersPerSecond];
                        %                     cnt2=cnt2+1;
                    end
                    
                end
            end
        end
    elseif length(line) > 1 & line(1:3) == 'Fix'
        Epoch = Epoch + 1;
        disp(['*** Epoch = ', num2str(Epoch)])
    end
    %% While 문 종료
    if line == -1, break;  end
end
fclose('all');