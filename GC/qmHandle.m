classdef qmHandle < handle
%   class qmHandle
%   DO : QM배열에서 원하는거 쉽게 빼오도록 만듬
%       생성자
%           qmHandle()
%       함수[함수명]                         [반환]            [설명]
%       1   getQM()                         QM              : 전체 QM 반환
%       2   pickQM(gs, prn, type)           QM              : 선택 QM 반환
%
%   Copyright: taeil Kim, July 21, 2015@INHA University
    
    properties
        QM;     % QM 배열
    end
    
    methods ( Access = public )
        %--- 생성자 -------------------------------------------------------
        function obj = qmHandle(QM)
        %                                                       @2015.07.21
            obj.QM = QM;
        end
        %--- 전체QM반환 ---------------------------------------------------
        function QM = getQM(obj)
        %                                                       @2015.07.21
            QM = obj.QM;
        end
        %--- 선택QM반환 ---------------------------------------------------
        function subQM = pickQM(obj, gs, prn, type)
        % input : gs    gps week second (범위연산 불가, 여러개 입력 가능)
        %         prn   PRN (범위연산 불가, 여러개 입력 가능)
        %         type  TYPE (범위연산 불가, 여러개 입력 가능)
        %                                                       @2015.07.21
            subQM = obj.QM;
            if gs ~= ':'    ,subQM = obj.pickRows(subQM, 1, gs);    end
            if prn ~= ':'   ,subQM = obj.pickRows(subQM, 2, prn);   end
            if type ~= ':'
                if type/100 < 1
                    subQM = obj.pickRows(subQM, 3, obj.allType(type));
                else
                    subQM = obj.pickRows(subQM, 3, type);
                end
            end
        end
    end
    
    methods ( Access = private )
        function table = pickRows(obj, table, column, operand)
            logic_ = zeros(length(table(:,1)), 1);
            for i = 1:length(operand)
                logic = logical ( table(:, column) == operand(i) );
                logic_ = logic_ + logic;
            end
            table  = table(logical(logic_), :);
        end
        function types = allType(obj, type)
            sys = [100; 200; 300];
            types = sys * logical(type) + logical(sys) * type;
            types = reshape(types, 1, []);
        end
    end
end