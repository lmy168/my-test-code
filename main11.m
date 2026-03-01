for i = 5.1
clearvars -except i;
goal=[30 30 30];              %目标点
source=[130 130 130];         %起始点  （source 来源）
GoalThreshold = 30;         % 设置目标点阈值
Delta =10;                  % 设置扩展步长
RadiusForNeib = 40;         % rewire的范围，半径r
MaxIterations = 2500;       % 最大迭代次数
UpdateTime = 50;            % 更新路径的时间间隔
DelayTime = 0.0;            % 绘图延迟时间
%% 建树初始化：T是树，v是节点
T.v(1).x = source(1);             % 把起始节点加入到T中
T.v(1).y = source(2); 
T.v(1).z = source(3);
T.v(1).xPrev = source(1);         % 节点的父节点坐标，起点的父节点是其本身
T.v(1).yPrev = source(2);
T.v(1).zPrev = source(3);
T.v(1).totalCost = 0;         % 从起始节点的累计cost，这里取欧式距离
T.v(1).indPrev = 0;           % 父节点的索引
% 绘制障碍物（以球为例，主要是方便计算）
x0 = 100; 
y0 = 100; 
z0 = 100;                         % 正方体中心坐标
cubeCenter = [95   130    15
    95   130    45
    95   130    75
    95   130   105
    95   130   135
    95   100    15
    95   100    45
    95   100    75
    95   100   105
    95   100    135
    65    130    15
    65    130    45
    65   100    15
    65   100    45
    100  40     15
    100  40     30
    130  40     15];       % 障碍物的中心
d = 30;                            % 正方体的边长为圆的直径

figure('color',[1 1 1]); % 创建白色背景的新图形窗口
figure(1); % 切换到图形窗口1
alphaValues = linspace(0, 0, 6); % 定义透明度渐变值
for k = 1:length(cubeCenter(:,1)) % i从1到5
    cubeVertices = d/2 * [-1 -1 -1; -1 -1 1; -1 1 1; -1 1 -1; 1 -1 -1; 1 -1 1; 1 1 1; 1 1 -1] + repmat(cubeCenter(k,:),8,1); % 计算每个正方体的顶点坐标
    cubeFaces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8]; % 定义每个正方体的面
    for j = 1:size(cubeFaces, 1)
        patch('Vertices', cubeVertices, 'Faces', cubeFaces(j,:), 'FaceColor', [1 1 1], 'FaceAlpha', alphaValues(j),'EdgeColor','none'); % 将蓝色改为橙色
    end
    hold on;
end
untitled3([80,115,0],[110,145,30])
untitled3([80,115,30],[110,145,60])
untitled3([80,115,60],[110,145,90])
untitled3([80,115,90],[110,145,120])
untitled3([80,115,120],[110,145,150])
untitled3([80,85,0],[110,115,30])
untitled3([80,85,30],[110,115,60])
untitled3([80,85,60],[110,115,90])
untitled3([80,85,90],[110,115,120])
untitled3([80,85,120],[110,115,150])
untitled3([50,115,0],[80,145,30])
untitled3([50,115,30],[80,145,60])
untitled3([50,85,0],[80,115,30])
untitled3([50,85,30],[80,115,60])
untitled3([85,25,0],[115,55,30])
untitled3([85,25,30],[115,55,60])
untitled3([115,25,0],[145,55,30])
axis equal; % 将当前图形的坐标轴比例设置为相等
hold on; % 保持当前的图形窗口
view(3); % 设置视角为三维
axis equal; % 将当前图形的坐标轴比例设置为相等
searchSize = [250 250 250]; % 搜索范围为250 250 250 搜索空间六面体
hold on; % 保持当前的图形窗口
view(3); % 设置视角为三维
%% 绘制起点和终点
hold on;                         
scatter3(source(1),source(2),source(3),'filled','g');  %scatter绘制三维散点图的函数，filled表示填充颜色为实心
scatter3(goal(1),goal(2),goal(3),'filled','b');        
tic;  % tic-toc: Functions for Elapsed Time            %开始计时，用于测量程序运行的时间
count = 1;                                   %将变量count初始化为1
pHandleList = [];                            %创建一个空的矩阵 pHandleList，用于存储图形绘制的句柄。
lHandleList = [];                            %创建一个空的矩阵 lHandleList，用于存储图形绘制的线句柄。
resHandleList = [];                          %创建一个空的矩阵 resHandleList，用于存储图形绘制的结果句柄。
findPath = 0;                                %将变量 findPath 初始化为 0，表示还未找到路径。
update_count = 0;                            %将变量 update_count 初始化为 0。
path.pos = [];                               %初始化path结构体中的pos字段为空矩阵，用于存储路径点的位置
% 绕着空间直线x=y=z，以十度为一个旋转分别得到每一个平面
k = [1, 1, 1]/sqrt(3); % 空间直线在xyz坐标系下的方向向量
jiaodu=i*10;
theta =jiaodu* pi / 180; % 旋转的角度，单位为弧度

% 定义需要旋转的平面方程
A = 0.13;
B = -0.76;
C = 0.63;
D = 0;

% 计算旋转矩阵
K = [0, -k(3), k(2); k(3), 0, -k(1); -k(2), k(1), 0]; % 叉乘矩阵
R = eye(3) + sin(theta) * K + (1 - cos(theta)) * K^2; % 罗德里格斯公式

% 将平面方程表示为向量形式
P = [A, B, C]';

% 计算旋转后的法向量
P_rotated = R' * P;

% 将旋转后的法向量转化为平面方程
P_rotated = P_rotated / norm(P_rotated);
D = 0;
A_rotated = P_rotated(1);
B_rotated = P_rotated(2);
C_rotated = P_rotated(3);
% 输出旋转前后的平面方程
fprintf('旋转前的平面方程：%.2fx + %.2fy + %.2fz + %.2f = 0\n', A, B, C, D);
fprintf('旋转后的平面方程：%.2fx + %.2fy + %.2fz + %.2f = 0\n', A_rotated, B_rotated, C_rotated, D);
A1 = A_rotated;
B1 = B_rotated;
C1 = C_rotated;
D = 0;
for iter = 1:MaxIterations                   %MaxIterations 是最大迭代次数，iter 是当前的迭代次数
    %Step 1:       在地图中随机采样一个点x_rand（sample）
if rand < 0.5
    while true
        % 生成一个三维坐标样本 x_rand
        x_rand = rand(1, 3).*searchSize;

        % 判断平面方程的值是否小于 0.01，如果是则跳出循环
        if abs(A1*x_rand(1) + B1*x_rand(2) + C1*x_rand(3) + D) < 0.01
            break;
        end
    end
else
    x_rand = goal;
    % 将目标点 goal 设置为样本 x_rand
    
end
    %Step 2:       遍历树，从树中找到最近邻近点x_near（Near）
    minDis = sqrt((x_rand(1) - T.v(1).x)^2 + (x_rand(2) - T.v(1).y)^2+(x_rand(3) - T.v(1).z)^2);
    minIndex = 1;               %给minlndex赋值为1，表示当前最小距离对应的索引为1
    for i = 2:size(T.v,2)	    %T.v按行向量存储，size（T.v，2）获得节点总数
    	distance = sqrt((x_rand(1) - T.v(i).x)^2 + (x_rand(2) - T.v(i).y)^2+(x_rand(3) - T.v(i).z)^2);   
        if(distance < minDis)
            minDis = distance;
            minIndex = i;   
        end     
    end
    x_near(1) = T.v(minIndex).x;    % 找到当前树中离x_rand最近的节点
    x_near(2) = T.v(minIndex).y;
    x_near(3) = T.v(minIndex).z;
    temp_parent = minIndex;         %将  minIndex 赋值给 temp_parent
    temp_cost = Delta + T.v(minIndex).totalCost;   %temp_cost是步长加上总代价
    %Step 3: 扩展得到x_new节点（Steer）
    movingVec = [x_rand(1)-x_near(1),x_rand(2)-x_near(2),x_rand(3)-x_near(3)];
    movingVec = movingVec/sqrt(sum(movingVec.^2));  % 单位化
    newPoint  = x_near + Delta * movingVec;         % 沿着当前最优方向移动一定距离得到新节点位置，以便在路径规划算法中生长路径树  （这里好像是矩阵的求法）
    %应该是得到新节点方便下面碰撞检测
    %plot3(x_rand(1), x_rand(2), x_rand(3),'ro', 'MarkerSize',10, 'MarkerFaceColor','r');
    %plot3(newPoint(1), newPoint(2), newPoint(3), 'bo', 'MarkerSize',10, 'MarkerFaceColor','b');
    
    %theta = atan2((x_rand(2) - x_near(2)),(x_rand(1) - x_near(1)));
    %x_new(1) = x_near(1) + cos(theta) * Delta;
    %x_new(2) = x_near(2) + sin(theta) * Delta;  
    %plot(x_rand(1), x_rand(2), 'ro', 'MarkerSize',10, 'MarkerFaceColor','r');
    %plot(x_new(1), x_new(2), 'bo', 'MarkerSize',10, 'MarkerFaceColor','b');
   
        
    % 检查节点是否是collision-free
    if ~checkPath3(x_near, newPoint, cubeCenter,d) % if extension of closest node in tree to the new point is feasible
        continue;
    end
 
    %Step 4: 在以x_new为圆心，半径为R的圆内搜索节点（NearC)    
    disToNewList = [];    % 每次循环要把队列清空
    nearIndexList = [];
    for index_near = 1:count            %count代表当前树节点数
        disTonew = sqrt((newPoint(1) - T.v(index_near).x)^2 + (newPoint(2) - T.v(index_near).y)^2+(newPoint(3) - T.v(index_near).z)^2);
        if(disTonew < RadiusForNeib)    % 满足条件：欧式距离小于R
            disToNewList  = [disToNewList disTonew];     % 满足条件的所有节点到x_new的cost
            nearIndexList = [nearIndexList index_near];     % 满足条件的所有节点基于树的索引
        end
    end
    
    %Step 5: 选择x_new的父节点，使x_new的累计cost最小（ChooseParent）  
    %这段代码的作用是在满足条件的节点中选择代价最小且路径可行的节点作为新节点的父节点，以构建路径树的生长。
    for cost_index = 1:length(nearIndexList)    % cost_index是基于disToNewList的索引，不是整棵树的索引
        costToNew = disToNewList(cost_index) + T.v(nearIndexList(cost_index)).totalCost;
        if(costToNew < temp_cost)    % temp_cost为通过minDist节点的路径的cost
            x_mincost(1) = T.v(nearIndexList(cost_index)).x;     % 符合剪枝条件节点的坐标
            x_mincost(2) = T.v(nearIndexList(cost_index)).y;
            x_mincost(3) = T.v(nearIndexList(cost_index)).z;
            if ~checkPath3(x_mincost, newPoint, cubeCenter,d) % if extension of closest node in tree to the new point is feasible
        continue;
            end
        	temp_cost = costToNew;
        	temp_parent = nearIndexList(cost_index);
        end
    end
    
    %Step 6: 将x_new插入树T（AddNodeEdge） 实现了路径树的构建和路径可视化功能
    count = count+1;    %  增加一个计数器，表示最新节点的索引
    
    T.v(count).x = newPoint(1);       %将最新节点的x，y，z坐标存储在树结构T的v字段中的第count个元素
    T.v(count).y = newPoint(2); 
    T.v(count).z = newPoint(3); 
    T.v(count).xPrev = T.v(temp_parent).x;  %将最新的父节点的x，y，z坐标存储在树结构T的v字段中的第count个元素的xPrev字段中
    T.v(count).yPrev = T.v(temp_parent).y;
    T.v(count).zPrev = T.v(temp_parent).z;
    T.v(count).totalCost = temp_cost;     %将最新节点的总代价存储在树结构T的v字段中的第count个元素的totalCost字段中
    T.v(count).indPrev = temp_parent;     %将最新节点的父节点的索引存储在树结构T的v字段中的第count个元素的indPrev字段中
    
   l_handle = plot3([T.v(count).xPrev; newPoint(1)], [T.v(count).yPrev; newPoint(2)],[T.v(count).zPrev; newPoint(3)], 'b', 'Linewidth', 1.3);  %在三维空间中绘制一条从父节点到最新节点的蓝色路径线段，并返回该路径线段的句柄
   p_handle = plot3(newPoint(1), newPoint(2),newPoint(3), 'ko', 'MarkerSize', 4, 'MarkerFaceColor','k');           %在三维空间中以黑色圆点形式绘制最新节点，并返回该节点的句柄
   
   pHandleList = [pHandleList p_handle];    %将最新节点的句柄添加到pHandleList列表中，用于后续的句柄管理和更新
   lHandleList = [lHandleList l_handle];    %将路径线段的句柄添加到lHandleList列表中，用于后续的句柄管理和更新
   pause(DelayTime);                %在绘图前暂停一段时间，以便观察路径的构建过程
    %Step 7: 剪枝(rewire)
    for rewire_index = 1:length(nearIndexList)            %对nearlndexList列表中的每个索引进行循环迭代，该列表包含与新节点（newPoint）相邻的节点的索引
        if(nearIndexList(rewire_index) ~= temp_parent)    % 若不是之前计算的最小cost的节点
            newCost = temp_cost + disToNewList(rewire_index);    % 计算新节点再返回路径的总代价。将temp_cost(从起点到父节点的总代价）与disToNewList（新节点到相邻节点的距离）相加
            if(newCost < T.v(nearIndexList(rewire_index)).totalCost)    % 检查新路径的总代价是否比相邻节点当前记录的总代价更小。是的话需要剪枝
                x_neib(1) = T.v(nearIndexList(rewire_index)).x;     % 符合剪枝条件节点的坐标
                x_neib(2) = T.v(nearIndexList(rewire_index)).y;
                x_neib(3) = T.v(nearIndexList(rewire_index)).z;
                if ~checkPath3(x_neib, newPoint, cubeCenter,d) % 检查从树中最接近新点的节点到新点的路径是否可行。如果不可行，则跳过当前迭代，继续下一个迭代
                    continue;
                end
                T.v(nearIndexList(rewire_index)).xPrev = newPoint(1);      % 对该neighbor信息进行更新
                T.v(nearIndexList(rewire_index)).yPrev = newPoint(2);
                T.v(nearIndexList(rewire_index)).zPrev = newPoint(3);
                T.v(nearIndexList(rewire_index)).totalCost = newCost;
                T.v(nearIndexList(rewire_index)).indPrev = count;       % 更新相邻节点的indPrev字段为新节点的索引
                delete(pHandleList());
                if nearIndexList(rewire_index) <= length(lHandleList)
                delete(lHandleList(nearIndexList(rewire_index)));
                end
                lHandleList(nearIndexList(rewire_index)) = plot3([T.v(nearIndexList(rewire_index)).x, newPoint(1)], [T.v(nearIndexList(rewire_index)).y, newPoint(2)],[T.v(nearIndexList(rewire_index)).z, newPoint(3)], 'r', 'Linewidth', 1.3);
                %在三维空间中绘制一条从相邻节点到新节点的红色路径线段，并更新lHandleList列表中相邻节点的路径线段句柄。
                %pHandleList = [pHandleList p_handle];    %绘图的句柄索引即为count
                %lHandleList = [lHandleList l_handle];
            end
        end
    end
    
    %Step 8:检查是否到达目标点附近
    disToGoal = sqrt((newPoint(1) - goal(1))^2 + (newPoint(2)- goal(2))^2+(newPoint(3)- goal(3))^2);
    %首先，代码计算新节点（newPoint）到目标点（goal）的距离（disToGoal）。如果距离小于目标阈值（GoalThreshold）且还没有找到路径，则进入if语句块。
    if(disToGoal < GoalThreshold && ~findPath)    % 找到目标点，此条件只进入一次
        findPath = 1;       %将findpath标记为1，表示已经找到了路径。
        count = count+1;    %增加计数器的值，用于标识新添加的目标节点在树中的索引。
        Goal_index = count;
        T.v(count).x = goal(1);          %将目标点的坐标goal赋值给树T中新节点的坐标
        T.v(count).y = goal(2);
        T.v(count).z = goal(3); 
        T.v(count).xPrev = newPoint(1);    %将新节点的父节点坐标赋值给树T中新节点的父节点坐标 
        T.v(count).yPrev = newPoint(2);
        T.v(count).yPrev = newPoint(3);
        T.v(count).totalCost = T.v(count - 1).totalCost + disToGoal;   %计算新节点到起始点的总代价，之前节点的总代价加上新节点到目标点的距离
        T.v(count).indPrev = count - 1;     %其父节点x_near的index     将新节点的父节点索引赋值为count-1
    end
    
    if(findPath == 1)                       %如果已经找到路径，则执行下面的代码块
        update_count = update_count + 1;    
        if(update_count == UpdateTime)      %将update_count等于Update Time,则表示需要更新路径
            update_count = 0;               
            j = 2;
            path.pos(1).x = goal(1);        %将目标点的坐标赋值给路径上的第一个位置（path.pos(1)）
            path.pos(1).y = goal(2);
            path.pos(1).z = goal(3);
            pathIndex = T.v(Goal_index).indPrev;    %通过Goal_index找到路径的起始点索引
            while 1                                 %通过while循环回溯路径，将路径的每个节点的坐标保存在path.pos中，直到回溯到起点（pathIndex == 0）为止。
                path.pos(j).x = T.v(pathIndex).x;
                path.pos(j).y = T.v(pathIndex).y;
                path.pos(j).z = T.v(pathIndex).z;
                pathIndex = T.v(pathIndex).indPrev;    % 沿终点回溯到起点
                if pathIndex == 0
                    break
                end
                j=j+1;
            end  
            
            for delete_index = 1:length(resHandleList)     %使用delete函数删除之前绘制的路径图形
            	delete(resHandleList(delete_index));
            end
            for j = 2:length(path.pos)                     %通过for循环遍历路径中的每个节点，绘制相应的路径线段，将图形句柄存储在resHandleList中
                res_handle = plot3([path.pos(j).x; path.pos(j-1).x;], [path.pos(j).y; path.pos(j-1).y],[path.pos(j).z; path.pos(j-1).z;], 'g', 'Linewidth', 4.5);
                resHandleList = [resHandleList res_handle];
            end
        end
    end  
	pause(DelayTime); %暂停DelayTime s，使得RRT*扩展过程容易观察
end
 
for delete_index = 1:length(resHandleList)       %使用for循环遍历resHandleList列表中的所有路径线段句柄，并调用delete函数将其删除
	delete(resHandleList(delete_index));
end
for j = 2:length(path.pos)           
	res_handle = plot3([path.pos(j).x; path.pos(j-1).x;], [path.pos(j).y; path.pos(j-1).y],[path.pos(j).z; path.pos(j-1).z;], 'g', 'Linewidth', 4.5);
	resHandleList = [resHandleList res_handle];
end
lastTotalCost = 0;  % 初始化变量用于存储满足条件的最后一个totalCost的值

for i = 1:length(T.v)
    if isnan(T.v(i).x)
        lastTotalCost = T.v(i).totalCost;  % 更新lastTotalCost为当前满足条件的totalCost的值
    end
end

fprintf('旋转角度：%d \n', jiaodu);
fprintf('满足条件的最后一个路径长度为：%f \n', lastTotalCost);     
disp('The path is found!');
elapsedTime = toc;   % 计算经过的时间
fprintf('经过的时间为：%f 秒\n\n', elapsedTime);
% 假设 pos 是包含 11 个结构体的 1x11 结构数组
pos = path.pos; % 替换为你的实际数据

% 获取结构数组的长度
num_pos = length(pos);

% 初始化矩阵，假设每个结构体中的字段名为 x 和 y
output_matrix = zeros(num_pos, 3);

% 展开结构数组的字段值
for p = 1:num_pos
    output_matrix(p, 1) = pos(p).x;
    output_matrix(p, 2) = pos(p).y;
    output_matrix(p, 3) = pos(p).z;
end

% 输出展开后的矩阵
 disp(output_matrix);
%%
% 定义矩阵
j = [output_matrix];

% 获取最后一列数据
lastColumn = j(:, end);

% 判断最后一列数据是否依次增大
if all(diff(lastColumn) >= 0)
    j = [output_matrix];
else
    % 找到第一个数据小于前面的数据的位置
    index = find(diff(lastColumn) < 0, 1);
    % 保留 z(z(:, end) < z(index+1, end), :) 的数据
    j1 = j(j(:, end) < j(index+1, end), :);
    % 保留 z(index+1:end, :) 的数据
    j2 = j(index+1:end, :);
    % 合并两个结果为一个矩阵
    j = [j1; j2];
end
 
 % 找到重复行的逻辑向量
duplicateRows = [false; all(diff(j) == 0, 2)];

% 删除重复行只保留第一次出现的行
j = j(~duplicateRows, :);
disp(j);
output_matrix=j ;
[j_rows, ~] = size(j);
disp(['矩阵 j 的行数为：' num2str(j_rows)]);
%% 增加距离的五次多项式
figure('color',[1 1 1]);
% 提供的数据
goal=[30 30 30];              %目标点
source=[130 130 130];         %起始点  （source 来源）
GoalThreshold = 30;         % 设置目标点阈值
Delta =10;                  % 设置扩展步长
RadiusForNeib = 40;         % rewire的范围，半径r
MaxIterations = 2500;       % 最大迭代次数
UpdateTime = 50;            % 更新路径的时间间隔
DelayTime = 0.0;            % 绘图延迟时间
%% 建树初始化：T是树，v是节点
T.v(1).x = source(1);             % 把起始节点加入到T中
T.v(1).y = source(2); 
T.v(1).z = source(3);
T.v(1).xPrev = source(1);         % 节点的父节点坐标，起点的父节点是其本身
T.v(1).yPrev = source(2);
T.v(1).zPrev = source(3);
T.v(1).totalCost = 0;         % 从起始节点的累计cost，这里取欧式距离
T.v(1).indPrev = 0;           % 父节点的索引
% 绘制障碍物（以球为例，主要是方便计算）
x0 = 100; 
y0 = 100; 
z0 = 100;                         % 正方体中心坐标
cubeCenter = [95   130    15
    95   130    45
    95   130    75
    95   130   105
    95   130   135
    95   100    15
    95   100    45
    95   100    75
    95   100   105
    95   100    135
    65    130    15
    65    130    45
    65   100    15
    65   100    45
    100  40     15
    100  40     30
    130  40     15]; 
d = 30;                            % 正方体的边长为圆的直径

% figure('color',[1 1 1]); % 创建白色背景的新图形窗口
% figure(2); % 切换到图形窗口1
alphaValues = linspace(0, 0, 6); % 定义透明度渐变值
for i = 1:length(cubeCenter(:,1)) % i从1到5
    cubeVertices = d/2 * [-1 -1 -1; -1 -1 1; -1 1 1; -1 1 -1; 1 -1 -1; 1 -1 1; 1 1 1; 1 1 -1] + repmat(cubeCenter(i,:),8,1); % 计算每个正方体的顶点坐标
    cubeFaces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8]; % 定义每个正方体的面
    for j = 1:size(cubeFaces, 1)
        patch('Vertices', cubeVertices, 'Faces', cubeFaces(j,:), 'FaceColor', [1 1 1], 'FaceAlpha', alphaValues(j),'EdgeColor','none'); % 将蓝色改为橙色
    end
    hold on;
end
untitled3([80,115,0],[110,145,30])
untitled3([80,115,30],[110,145,60])
untitled3([80,115,60],[110,145,90])
untitled3([80,115,90],[110,145,120])
untitled3([80,115,120],[110,145,150])
untitled3([80,85,0],[110,115,30])
untitled3([80,85,30],[110,115,60])
untitled3([80,85,60],[110,115,90])
untitled3([80,85,90],[110,115,120])
untitled3([80,85,120],[110,115,150])
untitled3([50,115,0],[80,145,30])
untitled3([50,115,30],[80,145,60])
untitled3([50,85,0],[80,115,30])
untitled3([50,85,30],[80,115,60])
untitled3([85,25,0],[115,55,30])
untitled3([85,25,30],[115,55,60])
untitled3([115,25,0],[145,55,30])
axis equal; % 将当前图形的坐标轴比例设置为相等
hold on; % 保持当前的图形窗口
view(3); % 设置视角为三维
axis equal; % 将当前图形的坐标轴比例设置为相等
searchSize = [250 250 250]; % 搜索范围为250 250 250 搜索空间六面体
hold on; % 保持当前的图形窗口
view(3); % 设置视角为三维

data =[output_matrix];

for i = 1:length(data(:,1))                       %对数据点遍历
    for j = 1:length(cubeCenter(:,1))           %对圆形障碍物进行循环遍历
        dist = norm(data(i,:) - cubeCenter(j,:)) - d/2;   %计算数据点与圆心之间的距离
        if dist < 15                                         %判断距离是否小于规定阈值（这里是15）
            newDist = d/2 + 14;                             % 增加一些距离
            dir = (data(i,:) - cubeCenter(j,:))/norm(data(i,:) - cubeCenter(j,:));
            data(i,:) = cubeCenter(j,:) + newDist*dir;
        end
    end
end
% for i = 1:length(data(:,1))                       %对数据点遍历
%     for j = 1:length(cubeCenter(:,1))           %对圆形障碍物进行循环遍历
%         dist = norm(data(i,:) - cubeCenter(j,:)) - sqrt((d^2) * 3);   %计算数据点与圆心之间的距离
%         if dist < 0                                         %判断距离是否小于规定阈值（这里是15）
%             newDist =sqrt((d^2) * 3)+ 3;                             % 增加一些距离
%             dir = (data(i,:) - cubeCenter(j,:))/norm(data(i,:) - cubeCenter(j,:));
%             data(i,:) = cubeCenter(j,:) + newDist*dir;
%         end
%     end
% end

X = data(:, 1);  % 获取自变量 X 的值
Y = data(:, 2);  % 获取第一个因变量 Y1 的值
Z = data(:, 3);  % 获取第二个因变量 Y2 的值

coefficients1 = polyfit(X, Y, 4);  % 对第一个因变量进行五次多项式拟合
x_range = linspace(X(1), X(end), 100); % 生成 x 的范围
fitted_curve1 = polyval(coefficients1, x_range);

% 绘制原始数据点
scatter3(X, Y, Z, 'filled', 'g');
hold on;

% 绘制拟合曲线
Z_fitted = (-A1*x_range - B1*fitted_curve1 ) / C1;

x_range_transposed = x_range';
fitted_curve1_transposed=fitted_curve1';
Z_fitted_transposed=Z_fitted';
% scatter3(x_range_transposed,fitted_curve1_transposed, Z_fitted_transposed, 'filled', 'b', 'SizeData', 20);
hold on;
% 生成一条三维线条
plot3(x_range_transposed,fitted_curve1_transposed, Z_fitted_transposed,  'r', 'LineWidth', 3);
hold on;

xlabel('X');
ylabel('Y');
zlabel('Z');
title('Fitted Polynomial Curves and 3D Line');
legend('Original Data', 'Fitted Curve', '3D Line', 'Location', 'best');
grid on;
hold on;
% 计算曲线长度
curve_length = 0; % 初始化曲线长度为0

for i = 1:length(x_range_transposed)-1
    % 计算相邻点之间的距离
    distance = sqrt((x_range_transposed(i+1) - x_range_transposed(i))^2 + (fitted_curve1_transposed(i+1) - fitted_curve1_transposed(i))^2 + (Z_fitted_transposed(i+1) - Z_fitted_transposed(i))^2);
    
    % 累加距离
    curve_length = curve_length + distance;
end
% 显示曲线长度
disp(['曲线长度为：', num2str(curve_length)]);
end
  