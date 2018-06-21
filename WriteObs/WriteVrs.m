function [] = WriteVrs(filename, yy, mo, dd)
%
%   input filename = VRS file
%         yy = year -> ex) 17
%         mo = month -> ex) 02
%         dd = day -> ex) 17
%
% coded by Joonseong Gim, mar 03, 2017
YY = yy + 2000;

fid = fopen(filename,'r');
if mo < 10
    fid_out = fopen(strcat('PTCO_',num2str(yy),'0',num2str(mo),num2str(dd),'_adm.txt'),'w');
else 
    fid_out = fopen(strcat('PTCO_',num2str(yy),num2str(mo),num2str(dd),'_adm.txt'),'w');
end

list = textscan(fid,'%s');
for i=1:length(list{1})
    line = list{1}{i};
    gga = strsplit(line,',');
    if length(gga) == 15
        FixQuality = str2num(cell2mat(gga(7)));
        %% UTC to gs
        utc = cell2mat(gga(2)); 
        h = str2num(utc(1:2));
        m = str2num(utc(3:4));
        s = str2num(utc(5:8));
        [gw gs] = date2gwgs(YY, mo, dd, h, m, s);
        %% GPGGA Latitude ddmm.mmmmmmmm to decimal degree
        Lati_ddmm = double(str2num(cell2mat(gga(3))));
        Lati_dd = fix(Lati_ddmm/100);
        Lati_mm = (Lati_ddmm*1e8 - Lati_dd*1e10)*1e-8 ;
        Latitude = Lati_dd + Lati_mm/60;
        %% GPGGA Longitude ddmm.mmmmmmmm to decimal degree
        Longi_ddmm = str2num(cell2mat(gga(5)));
        Longi_dd = fix(Longi_ddmm/100);
        Longi_mm = (Longi_ddmm*1e8 - Longi_dd*1e10)*1e-8 ;
        Longitude = Longi_dd + Longi_mm/60;
        %% GPGGA Height
        Height = str2num(cell2mat(gga(10))) + str2num(cell2mat(gga(12)));
        %% Decimal Degree to ECEF
        XYZ =  gd2xyz([Latitude, Longitude, Height]);
        %% 
        VRS(i,:) = [decimal(gs), Latitude, Longitude, Height, XYZ(1), XYZ(2), XYZ(3), FixQuality];
        fprintf(fid_out,'%6.2f %3.8d %3.8d %3.8d %7.8d %7.8d %7.8d \n',decimal(gs), Latitude, Longitude, Height, XYZ(1), XYZ(2), XYZ(3), FixQuality);
    end
    
end
