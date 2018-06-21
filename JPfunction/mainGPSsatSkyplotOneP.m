%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO: After reading rinex and Epoch, and then draw skyplot
% Input value
%   : file, epoch time(hhmmss)
%   : ex) mainGPSsatSkyplot('daej1350.12o', '132430')
% by Dong-Hyo Sohn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mainGPSsatSkyplot(O_file, epochTime)
clear all; clc; close all;
    O_file = 'rep21950.07o';
    epochTime = '185030';  % 23h 00m 00s
    
    site = strcat(O_file(1:4));
    doy = strcat(O_file(5:7));
    yy = strcat(O_file(10:11));
    
    % find navigation file
    Nname = strcat(site,doy,'0.',yy,'n'); % ex) 'daej1350.12n'
    N_file = ls(Nname);
    
    if exist(N_file) == 0
        Nname = strcat('brdc',doy,'0.',yy,'n');  % brdc file
        N_file = ls(Nname);
        
        if exist(N_file) == 0
            error('Navigation file not found.');
        end
    end
    % Make new file
    newfile = strcat(site,yy,doy,'-','AzEl.txt');
    
    delete(newfile);
    fid_Onew = fopen(newfile, 'a+');    
   
    %====================================================================
    %  [Read Navigation Message file and return an array named 'eph']
    eph = ReadEph(N_file);
    %--------------------------------------------------------------------

    %====================================================================
    %  [Read observation header part from observation file]
    fid_Obs = fopen(O_file, 'r');
    
    while 1
        str = fgets(fid_Obs);
        % [Get the approximate position from observation RINEX file]
        if str(61:67) == 'APPROX '              % The minimum number of char to a Head one-line is 67
            AppRcv = [str2num(str(1:14)) str2num(str(15:28)) str2num(str(29:42))];
            RcvX=AppRcv(1);
            RcvY=AppRcv(2);
            RcvZ=AppRcv(3);
        end
        % [Get number of observation types]
        if str(61:67) == '# / TYP'              % The minimum number of char to a Head one-line is 67
            nTyp = str2num(str(1:6));
            remTyp = ceil(nTyp/5);
        end
        % [End of Header]
        if str(61:67) == 'END OF '              % The minimum number of char to a Head one-line is 67
            break;
        end
    end
    %--------------------------------------------------------------------

    %====================================================================
    % [Calculate Latitute, Longitude, Height of the approximate position]
    % Output:  RcvLat - vector of ellipsoidal latitudes (radians)
    %          RcvLon - vector of ellipsoidal longitudes (radians)
    %          RcvHei - vector of ellipsoidal heights (m)
    [RcvLat, RcvLon, RcvHei] = xyz2ell(RcvX, RcvY, RcvZ);
    if isnan(RcvLat) || isnan(RcvLon) || isnan(RcvHei)
        error('Check [ APPROX POSITION XYZ ] of RINEX');
    end
    %--------------------------------------------------------------------
    
    %====================================================================
    % [Read body part from Observation file]
    newfile = strcat(site,yy,doy,'-','AzEl.txt');   %   Make new file
    
    
    delete(newfile);
    fid_Onew = fopen(newfile, 'a+');
    
    while ~feof(fid_Obs)
        % [Epoch by Epoch]
        str = fgets(fid_Obs);
        epochstr = str(1:26);
        newstr = str(1:29);

        [gw, gd, gs, doy] = calday2gpstime(str2num(str(2:3))+2000, str2num(str(5:6)), str2num(str(8:9)), str2num(str(11:12)), str2num(str(14:15)), str2num(str(16:26)));
        gts = gw*604800 + gs;
        ObsTime = gs;

        nSat = str2num(str(30:32));

        % Read time part of current epoch
        for i = 1 : nSat
            rsd = mod(i,12);
            cil = ceil(i/12);
            j = i - 12*(cil-1);

            if(cil > 1 && rsd == 1)
                str = fgets(fid_Obs);
            end

            satinfo(i,:) = str(3.*j+30:3.*j+32);    % save satellite PRN ('G','R','S','E')
            %sys(j) = str(3.*j+30);    % save navigation system ('G','R','S','E')
            epochPRN(i) = str2num(str(3.*j+31:3.*j+32));
        end

        % Read data part of current epoch
        AzEl = zeros(nSat, 2);
        for i = 1 : nSat            
            % Select only GPS data('G') in current epoch
            if(satinfo(i,1) == 'G')
                prn = str2num(satinfo(i,2:3));
                % Pick up the proper column in the ephemerides array
                icol = PickEph(eph, prn, gts);

                % Computes the XYZ coordinates of a given GPS satellite
                SatXYZ = GetSatPos(eph, icol, ObsTime);
                dPos(i,:) = SatXYZ - AppRcv;

                % Transform from XYZ Coord. to Topocentric Coord. using Latitude and Longitude
                topo(i,:) = xyz2topo(dPos(i,:), RcvLat, RcvLon);

                % Transform from Topocentric Coord. to Azimuth and Elevation
                AzEl(i,:) = topo2AzEl(topo(i,:));
                if AzEl(i,2) < 0
                    AzEl(i,2) = 0.0001;
                end				
            end
        end
        % Read data part
        for i = 1 : nSat
            flg = 1;
            typ = nTyp;
            for r = 1 : remTyp
                str = fgets(fid_Obs);                
                if(flg == 1 && typ < 6)
                    flg = 2;
                else
                    typ = typ - 5;
                end
            end
        end        

        % Write Epoch, PRN, Az, El, Signal strength
        for i = 1 : nSat
            fprintf(fid_Onew, '%26s %2d %8.4f %8.4f\n', epochstr, epochPRN(i), AzEl(i,1), AzEl(i,2));
            %fprintf(fid_Onew, '%26s %s%2d %8.4f %8.4f\n', epochstr, satinfo(i,1), epochPRN(i), AzEl(i,1), AzEl(i,2));
        end

        clear epochPRN;
        clear dPos;
        clear topo;
        clear AzEl;
        clear satinfo;
        clear nSat;
    end
    
    fclose(fid_Onew);
    fclose(fid_Obs);

    
%% ------------------------------------------------------------------------
% [ Plot one Epoch Sky Plot using Az, El ]
    SatAzEl = load(newfile);
    SatAzEl(:,10) = SatAzEl(:,4) + (SatAzEl(:,5)/60) + (SatAzEl(:,6)/(60*60));  % Time
    hr = str2num(epochTime(1:2));
    mi = str2num(epochTime(3:4));
    se = str2num(epochTime(5:6));
    epochT = unique(SatAzEl(:,10));
    
%--------------------------------------------------------------------------
%    Draw 1 epoch point Satellite position
    figure(1);
    hold on;    
%    Draw Circle
    circlea = 0:pi/30:2*pi;
    cx = cos(circlea);
    cy = sin(circlea);
%    Draw a major Circle 
    for i= [30 60 90]
        plot(cx*i, cy*i, '-', 'color', 'k', 'linewidth', 1);
    end
%    Draw a minor Circle
    for i=[15 45 75]
        plot(cx*i, cy*i, ':', 'color', 'k', 'linewidth', 1);
    end
%    Draw major Sector Lines inside a circle
    lenmax = 90;
    circleTick = (1:6)*pi/6;
    cosct = cos(circleTick); 
    sinct = sin(circleTick);
    cax = [-cosct; cosct];
    say = [-sinct; sinct];
    plot(lenmax*cax, lenmax*say, '-', 'color', 'k', 'linewidth', 1);
%    plot Satellite Az and El
    for i = 1 : length(epochT)
        if epochT(i) == (hr+(mi/60)+(se/(60*60)))            
            clear Az; clear El; clear xx; clear yy;
            
            itx = find(SatAzEl(:,10) == epochT(i));
            Az = SatAzEl(itx,8);
            El = SatAzEl(itx,9);
            xx = (El-90) .* -(sin(Az*pi/180));
            yy = (El-90) .* -(cos(Az*pi/180));
           
            plot(xx, yy, '.','color', 'r', 'Markersize', 20);
            text(xx+2, yy+2, num2str(SatAzEl(itx,7)), 'color', 'b', 'horizontalalignment', 'left', 'FontSize', 15);
            
            titNAME = strcat(newfile(1:9),'-',num2str(SatAzEl(itx(1),4),'%2d'),':',num2str(SatAzEl(itx(1),5),'%2d'),':',num2str(SatAzEl(itx(1),6),'%2d'));
            title(titNAME);
            break;
        end
    end
%    Insert direction text
    rlen = 1.06 * lenmax;
    for i = 1 : length(circleTick) 
        ticm1 = int2str(i*30);
        ticm2 = int2str(180+i*30);
         if ticm2 == '360'
             ticm2 =' ';
         end
        text( rlen*sinct(i),  rlen * cosct(i), ticm1, 'horizontalalignment', 'center');     
        text(-rlen*sinct(i), -rlen * cosct(i), ticm2, 'horizontalalignment', 'center');
    end
    set(gca,'FontWeight', 'bold');
    axis('equal');
    axis('off');
