close all

for i=1:7
    figure(i)
    switch i
        case 1
            
            tHour = mod(test1_average_result_vrs2(:,1), 86400);
            tHour = tHour/3600;
            
            subplot(3,5,[1,2,3,6,7,8,11,12,13])
            plot(test1_average_result_vrs2(:,3),test1_average_result_vrs2(:,2),'bo')
            grid on; hold on;
            plot([-5,5],[0,0],'r-')
            plot([0,0],[-5,5],'r-')
            axis([-5 5 -5 5])
            axis square
            xlabel({['H RMSE = ', num2str(rms(test1_average_result_vrs2(:,5)))],...
                ['3D RMSE = ',  num2str(rms(test1_average_result_vrs2(:,6)))]});
            
            subplot(3,5,[4,5])
            plot(tHour,test1_average_result_vrs2(:,2),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta N')
            xlabel('Hour')
            
            subplot(3,5,[9,10])
            plot(tHour,test1_average_result_vrs2(:,3),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta E')
            xlabel('Hour')
            
            subplot(3,5,[14,15])
            plot(tHour,test1_average_result_vrs2(:,4),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta V')
            xlabel('Hour')
            tHour = 0;
        case 2
            
             tHour = mod(test2_1_average_result_vrs2(:,1), 86400);
            tHour = tHour/3600;
            
            subplot(3,5,[1,2,3,6,7,8,11,12,13])
            plot(test2_1_average_result_vrs2(:,3),test2_1_average_result_vrs2(:,2),'bo')
            grid on; hold on;
            plot([-5,5],[0,0],'r-')
            plot([0,0],[-5,5],'r-')
            axis([-5 5 -5 5])
            axis square
            xlabel({['H RMSE = ', num2str(rms(test2_1_average_result_vrs2(:,5)))],...
                ['3D RMSE = ',  num2str(rms(test2_1_average_result_vrs2(:,6)))]});
            
            subplot(3,5,[4,5])
            plot(tHour,test2_1_average_result_vrs2(:,2),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta N')
            xlabel('Hour')
            
            subplot(3,5,[9,10])
            plot(tHour,test2_1_average_result_vrs2(:,3),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta E')
            xlabel('Hour')
            
            subplot(3,5,[14,15])
            plot(tHour,test2_1_average_result_vrs2(:,4),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta V')
            xlabel('Hour')
            tHour = 0;
        case 3
            
            tHour = mod(test2_2_average_result_vrs2(:,1), 86400);
            tHour = tHour/3600;
            
            subplot(3,5,[1,2,3,6,7,8,11,12,13])
            plot(test2_2_average_result_vrs2(:,3),test2_2_average_result_vrs2(:,2),'bo')
            grid on; hold on;
            plot([-5,5],[0,0],'r-')
            plot([0,0],[-5,5],'r-')
            axis([-5 5 -5 5])
            axis square
            xlabel({['H RMSE = ', num2str(rms(test2_2_average_result_vrs2(:,5)))],...
                ['3D RMSE = ',  num2str(rms(test2_2_average_result_vrs2(:,6)))]});
            
            subplot(3,5,[4,5])
            plot(tHour,test2_2_average_result_vrs2(:,2),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta N')
            xlabel('Hour')
            
            subplot(3,5,[9,10])
            plot(tHour,test2_2_average_result_vrs2(:,3),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta E')
            xlabel('Hour')
            
            subplot(3,5,[14,15])
            plot(tHour,test2_2_average_result_vrs2(:,4),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta V')
            xlabel('Hour')
            tHour = 0;
        case 4
            
            tHour = mod(test3_average_result_vrs2(:,1), 86400);
            tHour = tHour/3600;
            
            subplot(3,5,[1,2,3,6,7,8,11,12,13])
            plot(test3_average_result_vrs2(:,3),test3_average_result_vrs2(:,2),'bo')
            grid on; hold on;
            plot([-5,5],[0,0],'r-')
            plot([0,0],[-5,5],'r-')
            axis([-5 5 -5 5])
            axis square
            xlabel({['H RMSE = ', num2str(rms(test3_average_result_vrs2(:,5)))],...
                ['3D RMSE = ',  num2str(rms(test3_average_result_vrs2(:,6)))]});
            
            subplot(3,5,[4,5])
            plot(tHour,test3_average_result_vrs2(:,2),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta N')
            xlabel('Hour')
            
            subplot(3,5,[9,10])
            plot(tHour,test3_average_result_vrs2(:,3),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta E')
            xlabel('Hour')
            
            subplot(3,5,[14,15])
            plot(tHour,test3_average_result_vrs2(:,4),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta V')
            xlabel('Hour')
            tHour = 0;
        case 5            
            tHour = mod(test4_1_average_result_vrs2(:,1), 86400);
            tHour = tHour/3600;
            
            subplot(3,5,[1,2,3,6,7,8,11,12,13])
            plot(test4_1_average_result_vrs2(:,3),test4_1_average_result_vrs2(:,2),'bo')
            grid on; hold on;
            plot([-5,5],[0,0],'r-')
            plot([0,0],[-5,5],'r-')
            axis([-5 5 -5 5])
            axis square
            xlabel({['H RMSE = ', num2str(rms(test4_1_average_result_vrs2(:,5)))],...
                ['3D RMSE = ',  num2str(rms(test4_1_average_result_vrs2(:,6)))]});
            
            subplot(3,5,[4,5])
            plot(tHour,test4_1_average_result_vrs2(:,2),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta N')
            xlabel('Hour')
            
            subplot(3,5,[9,10])
            plot(tHour,test4_1_average_result_vrs2(:,3),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta E')
            xlabel('Hour')
            
            subplot(3,5,[14,15])
            plot(tHour,test4_1_average_result_vrs2(:,4),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta V')
            xlabel('Hour')
            tHour = 0;
        case 6
            tHour = mod(test4_2_average_result_vrs2(:,1), 86400);
            tHour = tHour/3600;
            
            subplot(3,5,[1,2,3,6,7,8,11,12,13])
            plot(test4_2_average_result_vrs2(:,3),test4_2_average_result_vrs2(:,2),'bo')
            grid on; hold on;
            plot([-5,5],[0,0],'r-')
            plot([0,0],[-5,5],'r-')
            axis([-5 5 -5 5])
            axis square
            xlabel({['H RMSE = ', num2str(rms(test4_2_average_result_vrs2(:,5)))],...
                ['3D RMSE = ',  num2str(rms(test4_2_average_result_vrs2(:,6)))]});
            
            subplot(3,5,[4,5])
            plot(tHour,test4_2_average_result_vrs2(:,2),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta N')
            xlabel('Hour')
            
            subplot(3,5,[9,10])
            plot(tHour,test4_2_average_result_vrs2(:,3),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta E')
            xlabel('Hour')
            
            subplot(3,5,[14,15])
            plot(tHour,test4_2_average_result_vrs2(:,4),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta V')
            xlabel('Hour')
            tHour = 0;
        case 7
            
            tHour = mod(test4_3_average_result_vrs2(:,1), 86400);
            tHour = tHour/3600;
            
            subplot(3,5,[1,2,3,6,7,8,11,12,13])
            plot(test4_3_average_result_vrs2(:,3),test4_3_average_result_vrs2(:,2),'bo')
            grid on; hold on;
            plot([-5,5],[0,0],'r-')
            plot([0,0],[-5,5],'r-')
            axis([-5 5 -5 5])
            axis square
            xlabel({['H RMSE = ', num2str(rms(test4_3_average_result_vrs2(:,5)))],...
                ['3D RMSE = ',  num2str(rms(test4_3_average_result_vrs2(:,6)))]});
            
            subplot(3,5,[4,5])
            plot(tHour,test4_3_average_result_vrs2(:,2),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta N')
            xlabel('Hour')
            
            subplot(3,5,[9,10])
            plot(tHour,test4_3_average_result_vrs2(:,3),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta E')
            xlabel('Hour')
            
            subplot(3,5,[14,15])
            plot(tHour,test4_3_average_result_vrs2(:,4),'b.:')
            grid on; hold on;
            xlim([tHour(1) tHour(length(tHour))]);
            ylim([-5 5])
            ylabel('\Delta V')
            xlabel('Hour')
            tHour = 0;
    end
end

% VRS ºñ±³ plot
for i=8:14
    figure(i)
    switch i
        case 8
            plot(test1_average_result(:,3),test1_average_result(:,2),'b.-')
            hold on; grid on;
            plot(test1_average_result_vrs(:,3),test1_average_result_vrs(:,2),'r.-')
            axis([-40 3 -21 8])
            axis equal
            legend('PPSoln','VRS')
        case 9
            plot(test2_1_average_result(:,3),test2_1_average_result(:,2),'b.-')
            hold on; grid on;
            plot(test2_1_average_result_vrs(:,3),test2_1_average_result_vrs(:,2),'r.-')
            axis([-40 3 -21 8])
            axis equal
            legend('PPSoln','VRS')
        case 10
            plot(test2_2_average_result(:,3),test2_2_average_result(:,2),'b.-')
            hold on; grid on;
            plot(test2_2_average_result_vrs(:,3),test2_2_average_result_vrs(:,2),'r.-')
            axis([-40 3 -21 8])
            axis equal
            legend('PPSoln','VRS')
        case 11
            plot(test3_average_result(:,3),test3_average_result(:,2),'b.-')
            hold on; grid on;
            plot(test3_average_result_vrs(:,3),test3_average_result_vrs(:,2),'r.-')
            axis([-60 -25 30 58])
            axis equal
            legend('PPSoln','VRS')
        case 12
            axis([-60 -25 30 58])
            plot(test4_1_average_result(:,3),test4_1_average_result(:,2),'b.-')
            axis([-60 -25 30 58])
            hold on; grid on;
            plot(test4_1_average_result_vrs(:,3),test4_1_average_result_vrs(:,2),'r.-')
            axis([-60 -25 30 58])
            axis equal
            axis([-60 -25 30 58])
            legend('PPSoln','VRS')
        case 13
            axis([-60 -25 30 58])
            plot(test4_2_average_result(:,3),test4_2_average_result(:,2),'b.-')
            axis([-60 -25 30 58])
            hold on; grid on;
            plot(test4_2_average_result_vrs(:,3),test4_2_average_result_vrs(:,2),'r.-')
            axis([-60 -25 30 58])
            axis equal
            axis([-60 -25 30 58])
            legend('PPSoln','VRS')
        case 14
            plot(test4_3_average_result(:,3),test4_3_average_result(:,2),'b-.')
            hold on; grid on;
            plot(test4_3_average_result_vrs(:,3),test4_3_average_result_vrs(:,2),'r.-')
            axis([-60 -25 30 58])
            axis equal
            legend('PPSoln','VRS')
    end
end

for i=1:14
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
        case 8
            figure(i)
            suptitle(['Sedan ', Method, '(Type A) ',Titlestr])
        case 9
            figure(i)
            suptitle(['Sedan ', Method, '(Type B-1) ',Titlestr])
        case 10
            figure(i)
            suptitle(['Sedan ', Method, '(Type B-2) ',Titlestr])
        case 11
            figure(i)
            suptitle(['SUV ', Method, '(Type A) ',Titlestr])
        case 12
            figure(i)
            suptitle(['SUV ', Method, '(Type B-1) ',Titlestr])
        case 13
            figure(i)
            suptitle(['SUV ', Method, '(Type B-2) ',Titlestr])
        case 14
            figure(i)
            suptitle(['SUV ', Method, '(Type B-3) ',Titlestr])
    end
end