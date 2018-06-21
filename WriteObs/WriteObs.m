function [] = WriteObs(obs_file)
%
%function [] = WriteObs(obs_files)
%
%   Read the given observation RINEX file and writes out 'QMfile'
%       Writing out the observables
%       GWS (STA) PRN OBS_TYPE OBS_VALUE
%
%   Originally coded by Jihye Won, May 3, 2006
%   Significantly modified by Kwan-Dong Park, July 23, 2006
%
% 1 = L1, 2 = L2, 3 = C1, 4 = P1, 5 = P2, 6 = D1, 7 = D2, 8 = S1, 9 = S2,
% 10 = C2,
% 111(L1), 112(L2), 120(C1), 121(P1), 122(P2), 123(C2), 131(D1), 132(D2), 141(S1),
% 142(S2)
%xj터
% Revised 9/6/13, KDPARK, 13개 이상의 위성을 감당하기 위해 일부 수정
% 2014.07.09 김미소
% Glonass도 Read할 수 있도록 수정 / GPS = TYPES, GLO = TYPES+100 현재는 GLO까지만 추가된 상태. 


%% 입출력 파일 정의
fid_obs = fopen(obs_file,'r');
fid_out = fopen('QMfile','w');
%% 관측치 종류 정의
TYPES = [111 112 120 121 122 131 132 141 142 123]';
[ObsSeq, NoTYPES] = GetObsSeq(obs_file);
Get2END(fid_obs);
%% ??
SatType = char('G','R','S','E'); % GPS = G, Glonass = R, SBAS = S, Galileo = 4
%% 관측치 출력 과정 --- c1(GS) c2(PRN) c3(OBS_TYPE) c4(OBS_MEASUREMENT)s
ready = 0;
while ready == 0
    s = fgets(fid_obs);
    if length(s) > 7
        Epoch(1,1:6) = sscanf(s,'%f',6)';
        if Epoch(1,1) > 80
            Y4 = Epoch(1,1) + 1900;
        else
            Y4 = Epoch(1,1) + 2000;
        end
        MM = Epoch(1,2); DDD = Epoch(1,3); HH = Epoch(1,4); MIN = Epoch(1,5); SS = Epoch(1,6);
        [gw,gs] = date2gwgs(Y4,MM,DDD,HH,MIN,SS);
        
        NoSats = str2num(s(31:32));
        arrPRN = zeros(NoSats,1);
        arrSatType = zeros(NoSats,1);
        
        fn = 34; fs = 33;
        for i = 1:NoSats
            if i == 13 | i == 25                 % Updated 9/6/2013, KDPARK, 13개 이상의 위성을 감당하기 위해
                s = fgets(fid_obs);     %
                fn = 34;                %
                fs = 33;
            end                         %
            PRN = str2num(s(fn:(fn + 1)));
            arrPRN(i) = PRN;
            fn = fn + 3;
            
            Type1 = s(fs);
            INDX_T1 = strmatch(Type1,SatType);
            arrSatType(i) = INDX_T1;
            fs = fs + 3;
        end
        
        if NoTYPES > 5
            for i = 1:NoSats
                s = fgets(fid_obs);
                sLen = length(s);
                ObsSeq1 = 1;
                fstr = 1;
                for j = 1:5
                    if (fstr + 13 < sLen)
                        Obs1 = str2num(s(fstr:(fstr + 13)));
                        if ~isempty(Obs1)
                            if arrSatType(i) == 1
                            fprintf(fid_out,'%8.1f %4d %5d %16.3f \n', gs, arrPRN(i), TYPES(ObsSeq(j)), Obs1);
                            else
%                             TYPES(ObsSeq(j)) = TYPES(ObsSeq(j)) + 100;
                            fprintf(fid_out,'%8.1f %4d %5d %16.3f \n', gs, arrPRN(i), TYPES(ObsSeq(j))+100, Obs1);
                            end
                        end
                        fstr = fstr + 16;
                        ObsSeq1 = ObsSeq1 + 1;
                    end
                end
                s = fgets(fid_obs);
                sLen = length(s);
                fstr = 1;
                for j = 6:NoTYPES
                    if (fstr + 13 < sLen)
                        Obs1 = str2num(s(fstr:(fstr + 13)));
                        if ~isempty(Obs1)
                             if arrSatType(i) == 1
                            fprintf(fid_out,'%8.3f %4d %5d %16.3f \n', gs, arrPRN(i), TYPES(ObsSeq(j)), Obs1);
                             else
%                              TYPES(ObsSeq(j)) = TYPES(ObsSeq(j)) + 100;
                            fprintf(fid_out,'%8.3f %4d %5d %16.3f \n', gs, arrPRN(i), TYPES(ObsSeq(j))+100, Obs1);
                            end
                        end
                        fstr = fstr + 16;
                        ObsSeq1 = ObsSeq1 + 1;
                    end
                end
            end
        else
            for i = 1:NoSats
                s = fgets(fid_obs);
                ObsSeq1 = 1;
                fstr = 1;
                sLen = length(s);
                for j = 1:NoTYPES
                    if fstr + 13 < sLen
                        Obs1 = str2num(s(fstr:(fstr + 13)));
                        if ~isempty(Obs1)
                            if arrSatType(i) == 1
                            fprintf(fid_out,'%8.3f %4d %5d %16.3f \n', gs, arrPRN(i), TYPES(ObsSeq(j)), Obs1);
                             else
%                              TYPES(ObsSeq(j)) = TYPES(ObsSeq(j)) + 100;
                            fprintf(fid_out,'%8.3f %4d %5d %16.3f \n', gs, arrPRN(i), TYPES(ObsSeq(j))+100, Obs1);
                            end
                        end
                        fstr = fstr + 16;
                        ObsSeq1 = ObsSeq1 +1;
                    end
                end
            end
        end
        
    else
        ready=1;
    end
end

fclose(fid_obs);
fclose(fid_out);
