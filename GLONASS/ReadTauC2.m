function [tau_c] = ReadTauC2(GLON_file)

% function [tau_c] = ReadTauC(GLON_file)
% GLONASS ��۱˵��� ���Ͽ��� TauC���� �о����
% 2014.10.12 ��̼�

% GLON_file = 'brdm2200.15p';
fid = fopen(GLON_file, 'r');

ready = 0;

while ~ready
    s = fgetl(fid);
    if length(s) <= 80
        if s(01 : 04) == 'GLUT'
            tau_c = str2num(s(06 : 22));
            break;
        end
    end
end