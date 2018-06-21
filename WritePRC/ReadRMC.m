function [year, month, day, hour, minute, second] = ReadRMC(s)
%
%function [year,month,day] = readRMC(string)
%
% DO: Extract Year-Month-Day from RMC
%
% <input>   string : RMC line
%
% <output>  Year(4-digit)/Month/Day : UTC
%
% <exception>
%           ReadRMC:NotEnoughData
%               ���� �Ǵٰ� �� ����(���� ���κ� �߻�)
%               <eg> $GPRMC,071509.0,A,3734.2632
%
% Copyright: Jinyi Kim, November 26, 2014@INHA University
%--- Modifications ---
% 11/27/2014 �����̰� ���� readRMC�� ���� ������
%  1/12/2015 taeil Kim, <NotEnoughData> ���� �߰�(1)
%  2/16/2015 taeil Kim, �ð��� �����ϵ��� ����

index=regexp(s(8:end),',','split'); % �޸������� ���ڿ� �и�
if length(index)~=12                % ������ ���� Ȯ��
    throw(MException('ReadRMC:NotEnoughData', ...
        'G#RMC:12���� �������� %d���� �����մϴ�.', length(index)));
end
%--- ��¥ ���� -------------------------------------------------------------
ddmmyy = index{9};
day   = str2num(ddmmyy(1:2));
month = str2num(ddmmyy(3:4));
year  = str2num(ddmmyy(5:end)) + 2000;
%--- �ð� ���� -------------------------------------------------------------
hhmmss = index{1};
hour  = str2num(hhmmss(1:2));
minute= str2num(hhmmss(3:4));
second= str2num(hhmmss(5:end));
end