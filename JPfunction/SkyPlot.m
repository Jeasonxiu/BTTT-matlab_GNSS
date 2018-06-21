classdef SkyPlot < handle
%   class SkyPlot
%   DO : SkyPlot 그리는 클래스, 자료처리 중에 데이터를 저장해서 속도향상
%       생성자
%           SkyPlot()
%       함수[함수명]                         [반환]            [설명]
%       1   add(prn, az, el, SNR)           -               : 데이터 삽입
%       2   PLOT()                     data [prn,az,el,SNR] : 그래프 그리기
%
%   Copyright: taeil Kim, July 27, 2015@INHA University

    properties
        SYSTEM = 'GBR';
        FIGURE = 277;
        iter   = 0;
        data;
        colorMap;
    end
    
    methods ( Access = public )
        function this = SkyPlot()
            this.data   = NaN(80000, 4);
        end
        
        function [] = add(this, prn, az, el, SNR)
            this.iter = this.iter + 1;
            this.data(this.iter, :) = ...
                [prn, fix(az*10^2)/10^2, fix(el*10^2)/10^2, SNR];
        end
        
        function data = PLOT(this)
            data = unique(this.data(1:this.iter, :), 'rows');
            this.skyplot_(data);
        end
    end
    
    methods ( Access = private )
        function [] = skyplot_(this, data)
            %--- SkyPlot의 가장 외곽원 그리기 ------------------------------
            figure(this.FIGURE);
            handle = newplot([]);
            hold(handle, 'on');
            this.colorMap = colormap(jet(60));
            rectangle('position', [-90, -90, 180, 180], ...
                'facecolor', 'white', 'Curvature', [1 1], 'linewidth', 2);
            %--- 방위각 범위 나누기 ----------------------------------------
            th = (1:6) * 2*pi / 12;
            cst = cos(th); snt = sin(th);
            cs = [cst; -cst];
            sn = [snt; -snt];
            %--- 방위각 표시 점선 그리기 (360/12) --------------------------
            line(90*sn, 90*cs, 'linestyle', ':', 'color', [0 0 0]);
            %--- 1.0 * 90은 최외곽원 위의 점이라 1.1 * 90에 텍스트 표기 -----
            rt = 1.1 * 90;
            %--- 방위각 텍스트 표기 ----------------------------------------
            for i = 1:length(th)
                text(rt*snt(i), rt*cst(i), int2str(i*30), ...
                    'horizontalalignment', 'center', 'FontWeight', 'bold');
                if i == length(th)
                    loc = num2str(0);
                else
                    loc = num2str(180 + i*30);
                end
                text(-rt*snt(i), -rt*cst(i), loc, ...
                    'horizontalalignment', 'center', 'FontWeight', 'bold');
            end
            
            %--- 고도각 표시 점선 원 그리기, 텍스트 표시하기 ----------------
            for elevation = 0 : 15 : 90
                elevationCircleRadius = 90*cos((pi/180) * elevation);
                
                rectangle('position', [-elevationCircleRadius, -elevationCircleRadius, ...
                    2*elevationCircleRadius, 2*elevationCircleRadius], ...
                    'lineStyle', ':', 'Curvature', [1 1]);
                
                text(0, elevationCircleRadius, num2str(elevation), ...
                    'BackgroundColor', 'white', 'horizontalalignment','center', 'FontWeight', 'bold');
            end
            axis([-95 95 -90 101]);
            %% 극좌표(고도각, 방위각)를 그래프 표시를 위해 직교좌표(X, Y)로 변환
            y = @(az, el) 90*cos(el * pi/180) .* cos(az * pi/180);
            x = @(az, el) 90*cos(el * pi/180) .* sin(az * pi/180);
            
            %% SkyPlot 그래프 그림 ================================
            %--- 위성의 위치 표시 ------------------------------------------
            for prn_ = unique(data(:,1))'
                data_ = data(data(:,1) == prn_, :);
                len   = length(data_(:,1));
                for m = 1:len
                    x_ = x(data_(m,2), data_(m,3));
                    y_ = y(data_(m,2), data_(m,3));
                    plot(handle, x_, y_, ...
                        '.', 'color', this.colorMap(data_(m,4), :));
                end
                plot(handle, x_, y_, 'o', 'color', this.colorMap(data_(len,4), :));
                plot(handle, x_, y_, '*', 'color', this.colorMap(data_(len,4), :));
                text(x_, y_, ...
                    ['  ' this.SYSTEM(fix(prn_/100)) num2str(mod(prn_,100))], ...
                    'color', 'b', 'FontWeight', 'bold');
            end
            colorbar;
            %--- 두 축의 길이 비율을 맞춤 ----------------------------------
            axis(handle, 'equal');
            %--- 기본 직교축 제거 ------------------------------------------
            axis(handle, 'off');
        end
    end
    
end

