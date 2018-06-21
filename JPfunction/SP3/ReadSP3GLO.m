function [arrSP3] = ReadSP3GLO(file_sp3)
%
% ReadSP3_GLO('infile.sp3')
% GLONASS�� sp3������ �о �ð��� ���� ��ķ� �����ϴ� �ڵ�
% result = [sec prn x(m) y(m) z(m) clk(10e-6s) yyyy mm dd hh min sec]
% 
% ������ ReadSP3�� GLONASS_SP3�� ����ǵ��� ������
% ����2 : ��̼�
%        GPS�� SP3���� -> GLONASS�� SP3�� �е��� ����(ReadSP3 -> ReadSP3_GLO)


%% �Լ� ����� �� input ����
% clc; clear all;
% file_sp3='igl18034.sp3';
% file_sp3 = 'Sta18034_sp3_glo.txt';
%
fid_sp3 = fopen(file_sp3,'r');
[s] = Get2ENDsp3(fid_sp3); %: ��� ���� -- �������� �ǵ����� ���?
%% ���ũ�� �̸� ���� �����ִ밹�� 24 X 96(15���� 24�ð��� 96��)
maxNsat = 24;
size_arrSP3 = maxNsat * 4;
arrSP3 = zeros(size_arrSP3, 6);
%% ���� ����
k = 0;
ready = 0;

while ~ready
    if s(1:3) == 'EOF'  %: ù ���� EOF�� ��ϵ� ��� ���� �����Ѵ�.
        break;
    end
    
    if s(1,1) == '*'    %: Time-Tagging���� "*"�� ����
        yr  = str2num(s(4:7));
        mon = str2num(s(9:10));
        day = str2num(s(12:13));
        hr  = str2num(s(15:16));
        min = str2num(s(18:19));
        sec = str2num(s(21:32));
        [gw, gs] = date2gwgs(yr, mon, day, hr, min, sec);
        
    elseif length(CutStr(s,' ',1)) == 2         %: PRN 1-9
        k = k + 1;
        prn = CutStr(s,' ',1);
        arrSP3(k,1) = round(gs);                %: time-tag (sec)
        arrSP3(k,2) = str2num(s(3:4));          %: PRN
        arrSP3(k,3) = str2num(s(6:18))*1000;    %: X (km->m)
        arrSP3(k,4) = str2num(s(19:32))*1000;   %: Y (km->m)
        arrSP3(k,5) = str2num(s(34:46))*1000;   %: Z (km->m)
        arrSP3(k,6) = str2num(s(47:60));        %: �ð���� (micro-sec)
        
    else length(CutStr(s,' ',1)) == 4;          %: PRN 10-24
        k = k + 1;
        prn = CutStr(s,' ',1);
        arrSP3(k,1) = round(gs);                %: time-tag (sec)
        arrSP3(k,2) = str2num(s(3:4));          %: PRN
        arrSP3(k,3) = str2num(s(6:18))*1000;    %: X (km->m)
        arrSP3(k,4) = str2num(s(19:32))*1000;   %: Y (km->m)
        arrSP3(k,5) = str2num(s(34:46))*1000;   %: Z (km->m)
        arrSP3(k,6) = str2num(s(47:60));        %: �ð���� (micro-sec)
        
    end
    
    s = fgets(fid_sp3);
    if feof(fid_sp3) == 1   %: ���� ���� �����ϸ� ����!
        ready = 1;
    end
end
fclose(fid_sp3);
