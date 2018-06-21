clear all
close all

%% VRS load
load('PTCO1_170116_adm.txt');
vrs= PTCO1_170116_adm;
vrs(:,1) = vrs(:,1)+18;

% loaddata = 'PTH_PP.mat';
% loaddata = 'PTH_PP_dop.mat';
% loaddata = 'PTH_Diff.mat';
loaddata = 'PTH_Diff_dop.mat';

load(loaddata);



TruePos = [-3058799.61420451,4083265.35912516,3814946.87192938];
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2);
cnt=1; 
for j=1:4
    if j == 1 || j == 3
        temp1 = (['test',num2str(j),'_L']);
        temp2 = (['test',num2str(j),'_R']);
        temp3 = (['test',num2str(j),'_average']);
    end
    for k=1:3
        if j==2
            if k < 3
                temp1 = (['test',num2str(j),'_',num2str(k),'_L']);
                temp2 = (['test',num2str(j),'_',num2str(k),'_R']);
                temp3 = (['test',num2str(j),'_',num2str(k),'_average']);
            end
        elseif j==4
            temp1 = (['test',num2str(j),'_',num2str(k),'_L']);
            temp2 = (['test',num2str(j),'_',num2str(k),'_R']);
            temp3 = (['test',num2str(j),'_',num2str(k),'_average']);
        end
        temp4 = (['FinalTTs=intersect(',temp1,'(:,1),',temp2,'(:,1));']);
        eval(temp4);
        for i=1:length(FinalTTs)
            gs = FinalTTs(i,1);
            vrs_coordi = vrs(find(vrs(:,1)==gs),5:7);
            
            temp5 =([temp3,'(i,1)= gs;']);
            temp6 =([temp3,'(i,2:4)=','(',temp1,'(find(',temp1,'(:,1)==gs),2:4)+',...
                temp2,'(find(',temp2,'(:,1)==gs),2:4))/2;']);
            temp11 = ([temp3,'(i,5:6)=[',temp1,'(i,7),',temp2,'(i,7)];']);
            eval(temp5);
            eval(temp6);
            eval(temp11);
            temp7 =(['GGA=',temp3,'(i,[2,3,4,5,6]);']);
            eval(temp7);
            dXYZ = [GGA(1), GGA(2), GGA(3)] - TruePos;
            
            Scount = [GGA(4:5)];
            dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
            
            dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
            
            dNE = sqrt(dN^2 + dE^2);        %rmsH = myRMS(dNE);
            
            d2D(i,1) = dNE;
            
            d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
            
            d3D(i,1) = d3;
            

            
            if ~isempty(vrs_coordi)
                
                dXYZ_vrs = vrs_coordi - TruePos;
                dXYZ_vrs2 = [GGA(1), GGA(2), GGA(3)] - vrs_coordi;
                dNEV_vrs = xyz2topo(dXYZ_vrs, TrueLat, TrueLon);
                gd2 = xyz2gd(vrs_coordi); TrueLat2 = gd2(1); TrueLon2 = gd2(2);
                dNEV_vrs2 = xyz2topo(dXYZ_vrs2, TrueLat2, TrueLon2);
                dN_vrs = dNEV_vrs(:,1); dE_vrs = dNEV_vrs(:,2); dV_vrs = dNEV_vrs(:,3);
                dN_vrs2 = dNEV_vrs2(:,1); dE_vrs2 = dNEV_vrs2(:,2); dV_vrs2 = dNEV_vrs2(:,3);
                dNE_vrs = sqrt(dN_vrs^2 + dE_vrs^2);        %rmsH = myRMS(dNE);
                dNE_vrs2 = sqrt(dN_vrs2^2 + dE_vrs2^2);        %rmsH = myRMS(dNE);
                d2D_vrs(i,1) = dNE_vrs;
                d2D_vrs2(cnt,1) = dNE_vrs2;
                d3_vrs = sqrt(dN_vrs.^2 + dE_vrs.^2 + dV_vrs.^2); %rms3 = myRMS(d3);
                d3_vrs2 = sqrt(dN_vrs2.^2 + dE_vrs2.^2 + dV_vrs2.^2); %rms3 = myRMS(d3);
                d3D_vrs(i,1) = d3_vrs;
                d3D_vrs2(cnt,1) = d3_vrs2;

                if dNE_vrs2 < 4
                    temp8 = ([temp3,'_result(i,:) = [gs, dN, dE, dV, dNE, d3, Scount];']);
                    eval(temp8)
                    temp13 = ([temp3,'_result_vrs(i,:) = [gs, dN_vrs, dE_vrs, dV_vrs, dNE_vrs, d3_vrs];']);
                    eval(temp13)
                    temp14 = ([temp3,'_result_vrs2(cnt,:) = [gs, dN_vrs2, dE_vrs2, dV_vrs2, dNE_vrs2, d3_vrs2];']);
                    eval(temp14)
                    cnt=cnt+1;
                end
            end
        end
        cnt=1;
    end
end

if loaddata(length(loaddata)-6:length(loaddata)-4) == 'dop'
    Titlestr = ' with Doppler Smoothing';
else
    Titlestr = ' w/o Doppler Smoothing';
end
if loaddata(7) == 'f'
    Method = 'Diff';
else
    Method = 'PP';
end



for i=1:7
    figure(i)
    switch i
        case 1
            
            tHour = mod(test1_average_result_vrs2(:,1), 86400);
            tHour = tHour/3600;
            subplot(4,4,[1,2,5,6])
            
            plot(test1_average_result(:,3),test1_average_result(:,2),'bo')
            hold on; grid on;
            plot(test1_average_result_vrs(:,3),test1_average_result_vrs(:,2),'ro')
            legend('PPSoln','VRS')
            axis ([-40 5 -30 15])
            axis square
            subplot(4,4,[9,10,13,14])
            plot(test1_average_result_vrs2(:,3),test1_average_result_vrs2(:,2),'bo')
            grid on; hold on;
            axis square
            xlabel({['H RMSE = ', num2str(rms(test1_average_result_vrs2(:,5)))],...
                ['3D RMSE = ',  num2str(rms(test1_average_result_vrs2(:,6)))]});
            subplot(4,4,[3,4])
            plot(tHour,test1_average_result_vrs2(:,2),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta N')
            xlabel('Hour')
            subplot(4,4,[7,8])
            plot(tHour,test1_average_result_vrs2(:,3),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta E')
            xlabel('Hour')
            subplot(4,4,[11,12])
            plot(tHour,test1_average_result_vrs2(:,4),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta V')
            xlabel('Hour')
            tHour = 0;
        case 2
            
            tHour = mod(test2_1_average_result_vrs2(:,1), 86400);
            tHour = tHour/3600;
            subplot(4,4,[1,2,5,6])
            
            plot(test2_1_average_result(:,3),test2_1_average_result(:,2),'bo')
            hold on; grid on;
            plot(test2_1_average_result_vrs(:,3),test2_1_average_result_vrs(:,2),'ro')
            legend('PPSoln','VRS')
            axis ([-40 5 -30 15])
            axis square
            subplot(4,4,[9,10,13,14])
            plot(test2_1_average_result_vrs2(:,3),test2_1_average_result_vrs2(:,2),'bo')
            grid on; hold on;
            axis square
            xlabel({['H RMSE = ', num2str(rms(test2_1_average_result_vrs2(:,5)))],...
                ['3D RMSE = ',  num2str(rms(test2_1_average_result_vrs2(:,6)))]});
            subplot(4,4,[3,4])
            plot(tHour,test2_1_average_result_vrs2(:,2),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta N')
            xlabel('Hour')
            subplot(4,4,[7,8])
            plot(tHour,test2_1_average_result_vrs2(:,3),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta E')
            xlabel('Hour')
            subplot(4,4,[11,12])
            plot(tHour,test2_1_average_result_vrs2(:,4),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta V')
            xlabel('Hour')
            tHour = 0;
        case 3
            
            tHour = mod(test2_2_average_result_vrs2(:,1), 86400);
            tHour = tHour/3600;
            subplot(4,4,[1,2,5,6])
            
            plot(test2_2_average_result(:,3),test2_2_average_result(:,2),'bo')
            hold on; grid on;
            plot(test2_2_average_result_vrs(:,3),test2_2_average_result_vrs(:,2),'ro')
            legend('PPSoln','VRS')
            axis ([-40 5 -30 15])
            axis square
            subplot(4,4,[9,10,13,14])
            plot(test2_2_average_result_vrs2(:,3),test2_2_average_result_vrs2(:,2),'bo')
            grid on; hold on;
            axis square
            xlabel({['H RMSE = ', num2str(rms(test2_2_average_result_vrs2(:,5)))],...
                ['3D RMSE = ',  num2str(rms(test2_2_average_result_vrs2(:,6)))]});
            subplot(4,4,[3,4])
            plot(tHour,test2_2_average_result_vrs2(:,2),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta N')
            xlabel('Hour')
            subplot(4,4,[7,8])
            plot(tHour,test2_2_average_result_vrs2(:,3),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta E')
            xlabel('Hour')
            subplot(4,4,[11,12])
            plot(tHour,test2_2_average_result_vrs2(:,4),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta V')
            xlabel('Hour')
            tHour = 0;
        case 4
            
            tHour = mod(test3_average_result_vrs2(:,1), 86400);
            tHour = tHour/3600;
            subplot(4,4,[1,2,5,6])
            
            plot(test3_average_result(:,3),test3_average_result(:,2),'bo')
            hold on; grid on;
            plot(test3_average_result_vrs(:,3),test3_average_result_vrs(:,2),'ro')
            legend('PPSoln','VRS')
            axis ([-65 -25 25 60])
            axis square
            subplot(4,4,[9,10,13,14])
            plot(test3_average_result_vrs2(:,3),test3_average_result_vrs2(:,2),'bo')
            grid on; hold on;
            axis square
            xlabel({['H RMSE = ', num2str(rms(test3_average_result_vrs2(:,5)))],...
                ['3D RMSE = ',  num2str(rms(test3_average_result_vrs2(:,6)))]});
            subplot(4,4,[3,4])
            plot(tHour,test3_average_result_vrs2(:,2),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta N')
            xlabel('Hour')
            subplot(4,4,[7,8])
            plot(tHour,test3_average_result_vrs2(:,3),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta E')
            xlabel('Hour')
            subplot(4,4,[11,12])
            plot(tHour,test3_average_result_vrs2(:,4),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta V')
            xlabel('Hour')
            tHour = 0;
        case 5
            
            tHour = mod(test4_1_average_result_vrs2(:,1), 86400);
            tHour = tHour/3600;
            subplot(4,4,[1,2,5,6])
            
            plot(test4_1_average_result(:,3),test4_1_average_result(:,2),'bo')
            hold on; grid on;
            plot(test4_1_average_result_vrs(:,3),test4_1_average_result_vrs(:,2),'ro')
            legend('PPSoln','VRS')
            axis ([-65 -25 25 60])
            axis square
            subplot(4,4,[9,10,13,14])
            plot(test4_1_average_result_vrs2(:,3),test4_1_average_result_vrs2(:,2),'bo')
            grid on; hold on;
            axis square
            xlabel({['H RMSE = ', num2str(rms(test4_1_average_result_vrs2(:,5)))],...
                ['3D RMSE = ',  num2str(rms(test4_1_average_result_vrs2(:,6)))]});
            subplot(4,4,[3,4])
            plot(tHour,test4_1_average_result_vrs2(:,2),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta N')
            xlabel('Hour')
            subplot(4,4,[7,8])
            plot(tHour,test4_1_average_result_vrs2(:,3),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta E')
            xlabel('Hour')
            subplot(4,4,[11,12])
            plot(tHour,test4_1_average_result_vrs2(:,4),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta V')
            xlabel('Hour')
            tHour = 0;
        case 6
            
            tHour = mod(test4_2_average_result_vrs2(:,1), 86400);
            tHour = tHour/3600;
            subplot(4,4,[1,2,5,6])
            
            plot(test4_2_average_result(:,3),test4_2_average_result(:,2),'bo')
            hold on; grid on;
            plot(test4_2_average_result_vrs(:,3),test4_2_average_result_vrs(:,2),'ro')
            legend('PPSoln','VRS')
            axis ([-65 -25 25 60])
            axis square
            subplot(4,4,[9,10,13,14])
            plot(test4_2_average_result_vrs2(:,3),test4_2_average_result_vrs2(:,2),'bo')
            grid on; hold on;
            axis square
            xlabel({['H RMSE = ', num2str(rms(test1_average_result_vrs2(:,5)))],...
                ['3D RMSE = ',  num2str(rms(test1_average_result_vrs2(:,6)))]});
            subplot(4,4,[3,4])
            plot(tHour,test4_2_average_result_vrs2(:,2),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta N')
            xlabel('Hour')
            subplot(4,4,[7,8])
            plot(tHour,test4_2_average_result_vrs2(:,3),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta E')
            xlabel('Hour')
            subplot(4,4,[11,12])
            plot(tHour,test4_2_average_result_vrs2(:,4),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta V')
            xlabel('Hour')
            tHour = 0;
        case 7
            
            tHour = mod(test4_3_average_result_vrs2(:,1), 86400);
            tHour = tHour/3600;
            subplot(4,4,[1,2,5,6])
            
            plot(test4_3_average_result(:,3),test4_3_average_result(:,2),'bo')
            hold on; grid on;
            plot(test4_3_average_result_vrs(:,3),test4_3_average_result_vrs(:,2),'ro')
            legend('PPSoln','VRS')
            axis ([-65 -25 25 60])
            axis square
            subplot(4,4,[9,10,13,14])
            plot(test4_3_average_result_vrs2(:,3),test4_3_average_result_vrs2(:,2),'bo')
            grid on; hold on;
            axis square
            xlabel({['H RMSE = ', num2str(rms(test4_3_average_result_vrs2(:,5)))],...
                ['3D RMSE = ',  num2str(rms(test4_3_average_result_vrs2(:,6)))]});
            subplot(4,4,[3,4])
            plot(tHour,test4_3_average_result_vrs2(:,2),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta N')
            xlabel('Hour')
            subplot(4,4,[7,8])
            plot(tHour,test4_3_average_result_vrs2(:,3),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta E')
            xlabel('Hour')
            subplot(4,4,[11,12])
            plot(tHour,test4_3_average_result_vrs2(:,4),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylabel('\Delta V')
            xlabel('Hour')
            tHour = 0;
    end
end
for i=1:7
    switch i
        case 1
            figure(i)
            suptitle(['Sedan ', Method, '(Type A) ',Titlestr])
        case 2
            figure(i)
            suptitle(['Sedan ', Method, '(Type B-1) ',Titlestr])
        case 3
            figure(i)
            suptitle(['Sedan ', Method, '(Type B-2) ',Titlestr])
        case 4
            figure(i)
            suptitle(['SUV ', Method, '(Type A) ',Titlestr])
        case 5
            figure(i)
            suptitle(['SUV ', Method, '(Type B-1) ',Titlestr])
        case 6
            figure(i)
            suptitle(['SUV ', Method, '(Type B-2) ',Titlestr])
        case 7
            figure(i)
            suptitle(['SUV ', Method, '(Type B-3) ',Titlestr])
    end
end