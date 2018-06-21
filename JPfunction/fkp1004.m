FileFKP = 'text_fkp_1004.txt';
in = load(FileFKP);
% c2(gw); c3(gs); c4(prn);
% c5(GPS L1 Code Indicator, 0 or 1)
% c6(GPS L1 Pseudorange) - RTCM 수신값 X 0.02
% c7(GPS L1 PhaseRange ? L1 Pseudorange) - RTCM 수신값 X 0.0005
% c8(GPS L1 Lock Time Indicator)
% c9(GPS Integer L1 Pseudorange Modulus Ambiguity) - RTCM 수신값 X 299792.458
% c10(GPS L1 CNR) - RTCM 수신값 X 0.25

%% 원시데이터 저장할 파일
fid_out = fopen('QYANP_15059','w');

Sats = unique(in(:,4));
NoSats = length(Sats);

tS = gs2h24(min(in(:,3)));
tE = gs2h24(max(in(:,3)));

for kS = 1:NoSats
    prn = Sats(kS);                 %: 특정위성 선택
    indxPRN = find(in(:,4) == prn);
    in1s = in(indxPRN, :);
    
    whichCol = 7;                   %: 그래프로 나타낼 칼럼 번호
    
    kE = length(in1s);
    out = zeros(kE, 4); %: 출력 c1...
    for k=1:kE
        out(k,1) = in1s(k, 3);
        out(k,2) = prn;
        out(k,3) = in1s(k, whichCol);
%         out(k,3) = in1s(k, 6) + in1s(k, 9);
        fprintf(fid_out,'%8d %3d 120 %12.3f\n', out(k,1), prn, out(k,3));
    end
    figure(prn)
    
    plot(gs2h24(out(:,1)), out(:,3), 'o:'); xlim([tS tE]);

end
%% 깔끔한 마무리
fclose(fid_out);