function [s] = Get2ENDsp3(fid_sp3)
%
%function [] = Get2ENDsp3(sp3_file)
% 
% DO: Get to the last part of SP3 Header
%
% <input> fid_sp3: SP3 File Handler
%
% Copyright: Kwan-Dong Park@LDEO October 25, 2013
%

%% SP3 ������ ����� #, +, %, /�� �Ǵ��ؼ� ����
ready = 0;
while ready == 0
    s = fgets(fid_sp3);
    switch s(1,1)
        case '#'
        case '+'
        case '%'
        case '/'
        otherwise
            ready = 1;  
    end
end
%% �� �� �ǵ����� ����� ��� ������� s�� ��� ������ ��
