classdef qmHandle < handle
%   class qmHandle
%   DO : QM�迭���� ���ϴ°� ���� �������� ����
%       ������
%           qmHandle()
%       �Լ�[�Լ���]                         [��ȯ]            [����]
%       1   getQM()                         QM              : ��ü QM ��ȯ
%       2   pickQM(gs, prn, type)           QM              : ���� QM ��ȯ
%
%   Copyright: taeil Kim, July 21, 2015@INHA University
    
    properties
        QM;     % QM �迭
    end
    
    methods ( Access = public )
        %--- ������ -------------------------------------------------------
        function obj = qmHandle(QM)
        %                                                       @2015.07.21
            obj.QM = QM;
        end
        %--- ��üQM��ȯ ---------------------------------------------------
        function QM = getQM(obj)
        %                                                       @2015.07.21
            QM = obj.QM;
        end
        %--- ����QM��ȯ ---------------------------------------------------
        function subQM = pickQM(obj, gs, prn, type)
        % input : gs    gps week second (�������� �Ұ�, ������ �Է� ����)
        %         prn   PRN (�������� �Ұ�, ������ �Է� ����)
        %         type  TYPE (�������� �Ұ�, ������ �Է� ����)
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