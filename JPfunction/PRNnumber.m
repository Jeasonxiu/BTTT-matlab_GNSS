function [SNN] = PRNnumber(S,NN)
%
% function [SNN] = PRNnumber(S,NN)
%
% DO: PRN Numbering (All GNSS system applied by RINEX V3.02)
%
% Copyright: Jinyi Kim, 17 Jan 2015, GPS Lab. INHA Univ.
%

% < 안내사항 >
% GPS, GLO, BDS 는 교수님 지정 넘버링
% GAL, QZS, SBAS 는 임의 넘버링 (추후 수정 필요)

switch S
    case 'G' % GPS                
        SN = 100;
    case 'C' % BEIDOU
        SN = 200;
    case 'R' % GLONASS
        SN = 300;
    case 'E' % GALILEO
        SN = 400;
    case 'J' % QZSS
        SN = 500;
    case 'S' % SBAS
        SN = 600;
    case 'I' % SBAS
        SN = 700;
end

SNN = SN+NN;