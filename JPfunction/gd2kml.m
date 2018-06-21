function gd2kml(filename, tt, gd, DOY,YY)

% SYNTAX:
%   gd2kml(filename, gs, gd, DOY,YY)
%
% INPUT:
%   filename = kml filename
%   gs = NX1 GPS week second
%   gd = NX3 Latitude, Longitude, Height
%   DOY = Day of year
%   YY = year(20xx)
%
% OUTPUT:
%   kml filename
%
% DESCRIPTION:
%   Write a KML file (Goole Earth) caculated position.

% clear all
% close all
% 
% filename = 'teheran';
% load('teheran_estm_2_sep_ja_ub.mat');
% tt = Jav_estm_gd(:,1);
% gd = Jav_estm_gd(:,2:4);
% DOY = 059;
% YY = 18;


%% 파일 생성
fkml=fopen([filename,'.kml'],'wt');

while (fkml == -1)
    fkml=fopen(kml_filename,'wt');
end
fprintf(fkml, '<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fkml, '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n');
fprintf(fkml, '<Document>\n');
fprintf(fkml, '\t<name>%s</name>\n', filename);
fprintf(fkml, '\t<open>1</open>\n');
fprintf(fkml, '\t<description>In Google Earth, enable historical imagery and sunlight, then click on each placemark to fly to that point in time.</description>\n');
fprintf(fkml, '\n');
fprintf(fkml, '\n');
fprintf(fkml, '\t<Folder>\n');
fprintf(fkml, '\t\t<name>%s Path with timestamps</name>\n',filename);
fprintf(fkml, '\t\t<open>1</open>\n');
fprintf(fkml, '\t\t<LookAt>\n');
fprintf(fkml, '\t\t\t<longitude>%s</longitude>\n',num2str(gd(1,2)));
fprintf(fkml, '\t\t\t<latitude>%s</latitude>\n',num2str(gd(1,1)));
fprintf(fkml, '\t\t\t<altitude>%s</altitude>\n',num2str(gd(1,3)));
fprintf(fkml, '\t\t\t<range>4060.093687469</range>\n');
fprintf(fkml, '\t\t\t<tilt>0</tilt>\n');
fprintf(fkml, '\t\t</LookAt>\n');
fprintf(fkml, '\t\t<Style id="pink-blank-dot">\n');
fprintf(fkml, '\t\t\t<IconStyle>\n');
fprintf(fkml, '\t\t\t\t<Icon>\n');
fprintf(fkml, '\t\t\t\t\t<href>http://maps.google.com/mapfiles/kml/paddle/pink-blank.png</href>\n');
fprintf(fkml, '\t\t\t\t</Icon>\n');
fprintf(fkml, '\t\t\t</IconStyle>\n');
fprintf(fkml, '\t\t</Style>\n');
%% 시작 - timestamp
for i = 1:length(tt)
    gs = tt(i);
    [gw, Gd] = ydoy2gwgd(YY, DOY);
    [yy, mo, dd, hh, mm, ss] = gwgs2date(gw, gs);
    if mo < 10
        mo = strcat('0',num2str(mo)); 
    else
        mo = num2str(mo);
    end
    if dd < 10
        dd = strcat('0',num2str(dd)); 
    else
        dd = num2str(dd);
    end
    if hh < 10
        hh = strcat('0',num2str(hh)); 
    else
        hh = num2str(hh);
    end
    if mm < 10
        mm = strcat('0',num2str(mm)); 
    else
        mm = num2str(mm);
    end
    ss= round(ss);
    if ss < 10
        ss = strcat('0',num2str(ss)); 
    else
        ss = num2str(ss);
    end
   
    longi = gd(i,2); lati = gd(i,1); alti = gd(i,3);
    Epoch = ['#',num2str(i)];
    
    
    fprintf(fkml, '\n');
    fprintf(fkml, '\t\t<Placemark>\n');
    fprintf(fkml, '\t\t<name>%s</name>\n',Epoch);
    fprintf(fkml, '\t\t\t<tessellate>1</tessellate>\n');
    fprintf(fkml, '\t\t\t<altitudeMode>absolute</altitudeMode>\n');
    fprintf(fkml, '\t\t\t<TimeStamp>\n');
    fprintf(fkml, '\t\t\t\t<when>%s-%s-%sT%s:%s:%s</when>\n',num2str(yy),mo,dd,hh,mm,ss);
    fprintf(fkml, '\t\t\t</TimeStamp>\n');
    fprintf(fkml, '\t\t\t<styleUrl>#pink-blank-dot</styleUrl>\n');
    fprintf(fkml, '\t\t\t<Point>\n');
    fprintf(fkml, '\t\t\t\t<coordinates>%s,%s,%s</coordinates>\n',num2str(longi),num2str(lati),num2str(alti));
    % fprintf(fkml, '\t\t\t\t<longitude>%s</longitude>\n',num2str(longi));
    % fprintf(fkml, '\t\t\t\t<latitude>%s</latitude>\n',num2str(lati));
    % fprintf(fkml, '\t\t\t\t<altitude>%s</altitude>\n',num2str(alti));
    fprintf(fkml, '\t\t\t</Point>\n');
    fprintf(fkml, '\t\t</Placemark>\n');
end
fprintf(fkml, '\t</Folder>\n');

%% 시작 - Path
fprintf(fkml, '\n');
fprintf(fkml, '\t<Placemark>\n');
fprintf(fkml, '\t\t<name>%s Path</name>\n',filename);
fprintf(fkml, '\t\t<Style>\n');
fprintf(fkml, '\t\t\t<LineStype>\n');
fprintf(fkml, '\t\t\t\t<color>ff0000ff</color>\n');
fprintf(fkml, '\t\t\t\t<width>2</width>\n');
fprintf(fkml, '\t\t\t</LineStype>\n');
fprintf(fkml, '\t\t</Style>\n');
fprintf(fkml, '\t\t<LineString>\n');
fprintf(fkml, '\t\t\t<tessellate>1</tessellate>\n');
% fprintf(fkml, '\t\t\t<altitudeMode>absolute</altitudeMode>\n');
fprintf(fkml, '\t\t\t<altitudeMode>clampToGround</altitudeMode>\n');
fprintf(fkml, '\t\t\t<coordinates>\n');
for i=1:length(tt)
    line = strcat(num2str(gd(i,2)),',',num2str(gd(i,1)),',',num2str(gd(i,3)));
%     line = strcat(num2str(gd(i,2)),',',num2str(gd(i,1)));
    fprintf(fkml, '\t\t\t%s\n',line);
end
fprintf(fkml, '\t\t\t</coordinates>\n');
fprintf(fkml, '\t\t</LineString>\n');
fprintf(fkml, '\t</Placemark>\n');

fprintf(fkml, '</Document>\n');
fprintf(fkml, '</kml>');