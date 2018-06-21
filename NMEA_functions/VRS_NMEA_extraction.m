function [] = VRS_NMEA_extraction(filename,yy,mo,dd);

% yy= 2016;mo=12;dd=02;
%% VRS_file load
vrsfile = filename;
% vrsfile = '1202T4.nmea';

%% VRS_NMEA_GPGGA ÃßÃâ

GPGGA_vrs = getGPGGA(vrsfile);
filename_vrs = vrsfile(1,1:max(length(vrsfile))-5);
filename_vrs = strcat(filename_vrs,'_GGA','.txt');
filename_vrs_fix = strcat(filename_vrs,'_GGA_fix','.txt');
filename_vrs_float = strcat(filename_vrs,'_GGA_float','.txt');
filename_vrs_dgps = strcat(filename_vrs,'_GGA_dgps','.txt');
filename_vrs_sa = strcat(filename_vrs,'_GGA_sa','.txt');
fileoutput_vrs = fopen(filename_vrs, 'w');
fileoutput_vrs_fix = fopen(filename_vrs_fix, 'w');
fileoutput_vrs_float = fopen(filename_vrs_float, 'w');
fileoutput_vrs_dgps = fopen(filename_vrs_dgps, 'w');
fileoutput_vrs_sa = fopen(filename_vrs_sa, 'w');

k = 1; ln = 1; fixline = 1; floatline = 1; saline = 1; dgpsline = 1;
for i = 1: length(GPGGA_vrs)
    vrs_GGA = GPGGA_vrs{i,1};
    
    line = char(vrs_GGA);
    fprintf(fileoutput_vrs, '%s\r\n', line);
    if length(vrs_GGA) > 8
        [hh,mm,ss,x,y,z,la,lo,qi,nSats,ht] = NEWreadGGA3(vrs_GGA) ;
        [gw, gs] = date2gwgs(yy, mo, dd, hh, mm, ss);
        gpgga_vrs(k,:) = [round(gs),x,y,z,la,lo,qi,nSats,ht,ln];
        k = k + 1;
    end
    if qi == 4
        Fixline(fixline,:) = [round(gs),ln];
        fixline = fixline + 1;
        fprintf(fileoutput_vrs_fix, '%s\r\n', line);
    elseif qi == 5
        Floatline(floatline,:) = [round(gs),ln];
        floatline = floatline + 1;
        fprintf(fileoutput_vrs_float, '%s\r\n', line);
    elseif qi == 1
        SAline(saline,:) = [round(gs),ln];
        saline = saline + 1;
        fprintf(fileoutput_vrs_sa, '%s\r\n', line);
    elseif qi == 2
        DGPSline(dgpsline,:) = [round(gs),ln];
        dgpsline = dgpsline + 1;
        fprintf(fileoutput_vrs_dgps, '%s\r\n', line);
    end
    ln = ln + 1;
end