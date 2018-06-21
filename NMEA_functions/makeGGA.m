% function [] = makeGGA(estm, year, month, day)
%
%   [] = makeGGA(estm, year, month, day)
%   
%   input estm  : Position result(n X 5 or n X more)
%               ex> [gs x y z .....]
%         year  : UTC������ ���� �⵵
%         month : UTC������ ���� ��
%         day   : UTC������ ���� ��
%
%   output nmea text file
%               ex> GPGGA_yearmonthday.txt
%
%   example> makeGGA(estm,2017,03,28);



clear all
close all
load('DD_170216_kf.mat')
year =2017; month=03;day=28;
estm = Base;

%% �Է� ����, ��, ���� �߸� �Ǿ��� ��� ����
stop = 0;
if month > 13
    stop = 1;
else
    if month == 1 | month == 3 | month == 5 | month == 7 | month == 8 | month == 10 | month == 12
        if day >= 32
            stop = 2;
        end
    elseif month == 4 | month == 6 | month == 9 | month == 11
        if day >= 31
            stop = 2;
        end
    elseif month == 2
        if year == 2004 | year == 2008 | year == 2012 | year == 2016 | year == 2020 | year == 2024
            if day >= 30
                stop = 2;
            end
        else
            if day >= 29
                stop = 2;
            end
        end
    end
end
if stop == 1
    disp('�� �Է��� �߸��Ǿ����ϴ�.')
elseif stop == 2
    disp('�� �Է��� �߸��Ǿ����ϴ�.')
end
if stop == 0
    filename = strcat('GPGGA_',num2str(year),num2str(month),num2str(day),'.txt');
    fid_out = fopen(filename,'w');
    [gw gs]=date2gwgs(year, month, day, 0, 0, 0);       % date�� ������� gw ����
    
    for i = 1:length(estm)
        % for i = 10:10
        gs = estm(i,1) - 18;
        xyz = estm(i,2:4);
        %% UTC�� ����� ����
        [yy, mo, dd, hh, mm, ss] = gwgs2date(gw, gs);       % gw�� gs�� ������ date ����
        hh = round(hh)-9; mm = round(mm); ss = round(ss);
        if ss < 10
            ss = strcat('0',num2str(ss));
        elseif ss == 60             % 60�ʰ� �Ǹ� ���� �ø�
            mm = mm + 1;
            ss = '00';
        else
            ss = num2str(ss);
        end
        if mm < 10
            mm = strcat('0',num2str(mm));
        elseif mm == 60             % 60���� �Ǹ� �ð��� �ø�
            hh = hh + 1;
            mm = '00'
        else
            mm = num2str(mm);
        end
        if hh < 10
            hh = strcat('0',num2str(hh));
        elseif hh == 24             % 24�ð��� �Ѿ�� �������� �ø�
            dd = dd + 1;
            hh == '00'
        else
            hh = num2str(hh);
        end
        UTC = strcat(hh, mm, ss,'.00');
        gd = xyz2gd(xyz);
        Lati_ddmm = fix(gd(1))*100 + ((gd(1)-fix(gd(1)))*60);
        Longi_ddmm = fix(gd(2))*100 + ((gd(2)-fix(gd(2)))*60);
        Altitude = gd(3);
        %% GPGGA ������ ���� ����
        header = '$GPGGA'; North = 'N'; East = 'E';Fix_quality = '1'; NumSats = '10';
        geoid = '00.0';HDOP = '0.00';Meter = 'M';tale = '*77';
        %% GPGGA �ؽ�Ʈ�� ����� ����
        fprintf(fid_out,'%s,%s,%4.5f,%s,%5.5f,%s,%s,%s,%s,%3.1f,%s,%s,%s,,%s \n',header, UTC,...
            Lati_ddmm, North, Longi_ddmm, East, Fix_quality, NumSats, HDOP, Altitude, Meter, geoid, Meter, tale);
    end
end