function GGA = RTAP2GGA(filename)
%
%   function GGA = RTAP2GGA(filename)
%
%   input filename = rtap.txt filename
%
%   output GGA : NX1 cell array & GGA.txt
%
% coded by Joonseong GIM, Feb, 14, 2017

% clear all
% 
% filename = 'T3_RTAP_161207.txt';

%% filename load
data = load(filename);

%% GPGGA.txt »ý¼º
fid_out = fopen('GPGGA.txt', 'w');
for i = 1:length(data(:,1))
    %% ex> $GPGGA,042603.00,3739.63426490,N,12644.98709341,E,4,15,0.7,12.93343,M,18.79776,M,1.0,0859*75
    header = '$GPGGA,';
    utc = num2str(data(i,1));
    if length(utc) <= 5
        UTC = strcat('0',num2str(data(i,1)),'.00,');
    else
        UTC = strcat(num2str(data(i,1)),'.00,');
    end
    
    Lati_ddmm = num2str(deg2dm(data(i,2)));
    if length(Lati_ddmm) == 7
        Lati_ddmm = strcat(num2str(deg2dm(data(i,2))),'00');
    elseif length(Lati_ddmm) == 8
        Lati_ddmm = strcat(num2str(deg2dm(data(i,2))),'0');
    end
    Longi_ddmm = num2str(deg2dm(data(i,3)));
    if length(Longi_ddmm) == 8
        Longi_ddmm = strcat(num2str(deg2dm(data(i,3))),'00');
    elseif length(Longi_ddmm) == 9
        Longi_ddmm = strcat(num2str(deg2dm(data(i,3))),'0');
    end
    north = ',N,'; east = ',E,';
    height = num2str(data(i,4)-1.19);
    if length(height) == 4
        height = strcat(num2str(data(i,4)-1.19),'0000');
    elseif length(height) == 5
        height = strcat(num2str(data(i,4)-1.19),'000');
    elseif length(height) == 6
        height = strcat(num2str(data(i,4)-1.19),'00');
    elseif length(height) == 2
        height = strcat(num2str(data(i,4)-1.19),'.00000');
    end
   
    remains = strcat('4,15,0.7,12.93343,M,',height,',M,1.19,0859*75');
    line = strcat(header, UTC, Lati_ddmm, north, Longi_ddmm, east, remains);
    fprintf(fid_out, '%s \r\n',line);
end


