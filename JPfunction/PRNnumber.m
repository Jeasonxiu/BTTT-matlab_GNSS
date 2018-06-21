function [SNN] = PRNnumber(S,NN)
%
% function [SNN] = PRNnumber(S,NN)
%
% DO: PRN Numbering (All GNSS system applied by RINEX V3.02)
%
% Copyright: Jinyi Kim, 17 Jan 2015, GPS Lab. INHA Univ.
%

% < �ȳ����� >
% GPS, GLO, BDS �� ������ ���� �ѹ���
% GAL, QZS, SBAS �� ���� �ѹ��� (���� ���� �ʿ�)

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