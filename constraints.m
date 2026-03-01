function [c, ceq] = constraints(C)
    c = [1 - reshape(C, 4, [])]; % 约束条件：C 需要大于等于 1
    ceq = [];
end
%总得来说，这段代码定义了一个简单的约束条件，要求向量C中的所有元素都大于等于1。