% function feasible=feasiblePoint3(point,circleCenter,r) %定义了一个函数feasiblePoint3，point待检查的点的坐标
%                             %circleCenter包含多个圆心坐标的矩阵，每行代表一个圆的圆心坐标；
%                             %r包含多个圆半径的向量，每个元素对应一个圆的半径
% feasible=true;              %设置变量feasible的初始值为true，表示默认情况下认为给定点是可行的
% % check if collission-free spot and inside maps
%     for row = 1:length(circleCenter(:,1))      %对每个圆的圆心进行遍历
%         if sqrt(sum((circleCenter(row,:)-point).^2)) <= r(row)   %计算当前圆心和给定点之间的欧氏距离，并将其与对应圆的半径进行比较。
%         feasible = false;                      %如果给定点不在任何圆内，则保持feasible为true，表示给定点是安全的。
%         break;
%         end
%     end
% end

function feasible = feasiblePoint3(point, cubeCenter, d)
    feasible = true; % 默认情况下认为给定点是可行的
    
    % 检查给定点是否在任何一个正方体内部
    for row=1:length(cubeCenter(:,1))
        cubeMin = cubeCenter(row,:) - d/2; % 计算正方体的最小坐标
        cubeMax = cubeCenter(row,:) + d/2; % 计算正方体的最大坐标
        
        % 判断给定点是否在当前正方体内部
        if point(1) >= cubeMin(1) && point(1) <= cubeMax(1) && ...
           point(2) >= cubeMin(2) && point(2) <= cubeMax(2) && ...
           point(3) >= cubeMin(3) && point(3) <= cubeMax(3)
            feasible = false;
            break;
        end
    end
end
