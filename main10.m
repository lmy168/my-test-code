%% 增加距离的五次多项式
figure('color',[1 1 1]);
% 提供的数据
goal=[10 10 10];              %目标点
source=[150 150 150];         %起始点  （source 来源）
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
cubeCenter = [95   100    15
    95   100    45
    95   100    75
    95   100   105
    70    50    15
    70    50    45
    70    50    75
    65   100    15
    65   100    45];    % 障碍物的中心
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
untitled3([80,85,0],[110,115,30])
untitled3([80,85,30],[110,115,60])
untitled3([80,85,60],[110,115,90])
untitled3([80,85,90],[110,115,120])
untitled3([55,35,0],[85,65,30])
untitled3([55,35,30],[85,65,60])
untitled3([55,35,60],[85,65,90])
untitled3([50,85,0],[80,115,30])
untitled3([50,85,30],[80,115,60])
axis equal; % 将当前图形的坐标轴比例设置为相等
hold on; % 保持当前的图形窗口
view(3); % 设置视角为三维
axis equal; % 将当前图形的坐标轴比例设置为相等
searchSize = [250 250 250]; % 搜索范围为250 250 250 搜索空间六面体
hold on; % 保持当前的图形窗口
view(3); % 设置视角为三维

data =[  10.0000   10.0000   10.0000
   13.9172   13.8576   10.3144
   28.8472   29.0825   43.3097
   43.9441   44.3713   71.1936
   59.1703   59.8704  101.3406
   75.1585   75.7807  114.1079
   95.5880   96.1501  130.7873
  119.5621  119.8553  137.9564
  131.6717  131.8145  140.6565
  150.0000  150.0000  150.0000];

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