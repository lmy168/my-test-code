%% checkPath3.m	         %定义一个函数checkPath3
function feasible=checkPath3(n,newPos,circleCenter,r)        %该函数接受四个输入参数，分别是当前节点位置n
                                                             %新位置newPos、圆心位置circleCenter和半径r
feasible=true;              %初始化变量feasible为true，假设路径是可行的
movingVec = [newPos(1)-n(1),newPos(2)-n(2),newPos(3)-n(3)];   %计算从当前节点到新位置的移动向量，即新位置减去当前位置得到的差值
movingVec = movingVec/sqrt(sum(movingVec.^2));        %将移动向量进行单位化，除以向量的模长，使得向量成为单位向量，方便后续使用
for R=0:0.5:sqrt(sum((n-newPos).^2))                  %遍历从当前节点到新位置之间的距离范围，步长为0.5，范围是从0到当前节点到新位置的欧氏距离
    posCheck=n + R .* movingVec;         %根据当前节点、移动向量和距离R，计算出当前检查点的位置。
    if ~(feasiblePoint3(ceil(posCheck),circleCenter,r) && feasiblePoint3(floor(posCheck),circleCenter,r))  %调用函数feasiblePoint3来检查当前检查点的上方和下方是否可行。如果上方和下方至少有一个不可行，则执行下面的代码。
        feasible=false;break;            %将变量feasible设置为false，表示路径不可行，并跳出循环。
    end
end
    if ~feasiblePoint3(newPos,circleCenter,r), feasible=false;     %检查新位置是否可行，如果不可行，将变量feasible设置为false。
    end
end