function [tau_c] = ReadTauC(GLON_file)

% function [tau_c] = ReadTauC(GLON_file)
% GLONASS ��۱˵��� ���Ͽ��� TauC���� �о����
% 2014.10.12 ��̼�

% GLON_file = 'brdc2760.14g';
fid = fopen(GLON_file, 'r');

ready = 0;

while ~ready
    s = fgetl(fid);
    if length(s) <= 80
        if s(61 : 64) == 'CORR'
            tau_c = str2num(s(23 : 40));
            break;
        end
    end
end