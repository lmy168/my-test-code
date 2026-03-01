% 清理并初始化
clearvars;
close all;
clc;

% 起始点和目标点
source = [150 150 150];
goal = [10 10 10];

% 平面方程系数
A = 0.13;
B = -0.76;
C = 0.63;
D = 0;

% 旋转轴（x=y=z直线）
k = [1, 1, 1] / sqrt(3);

% 清理并初始化
clearvars;
close all;
clc;

% 起始点和目标点
source = [150 150 150];
goal = [10 10 10];

% 平面方程系数
A = 0.13;
B = -0.76;
C = 0.63;
D = 0;

% 旋转轴（x=y=z直线）
k = [1, 1, 1] / sqrt(3);

% 创建图形窗口
figure('color',[1 1 1]);
hold on;
grid on;
axis equal;
view(3);
xlabel('X');
ylabel('Y');
zlabel('Z');
%title('起始点、终点、轴线与旋转平面');

% 绘制起始点 - 使用柔和的绿色
scatter3(source(1), source(2), source(3), 50, 'filled', ...
    'MarkerFaceColor', [0.3 0.7 0.3], 'DisplayName', '起始点'); % 柔和的绿色

% 绘制终点 - 保持红色
scatter3(goal(1), goal(2), goal(3), 50, 'filled', ...
    'MarkerFaceColor', [0.8 0.2 0.2], 'DisplayName', '终点'); % 柔和的红色

% 绘制连接起始点和终点的轴线 - 使用灰蓝色
plot3([source(1), goal(1)], [source(2), goal(2)], [source(3), goal(3)], ...
      'b-', 'LineWidth', 2, 'Color', [0.8 0.2 0.2], 'DisplayName', '轴线');

% 定义不同的颜色（9种不同颜色）
colors = [
    0.6 0.8 1;    % 浅蓝色
    0.8 0.6 1;    % 浅紫色
    0.6 1.0 0.8;  % 浅绿色
    1.0 0.8 0.6;  % 浅橙色
    0.8 1.0 0.6;  % 黄绿色
    1.0 0.6 0.8;  % 浅粉色
    0.6 0.9 1.0;  % 天蓝色
    0.9 0.6 1.0;  % 紫粉色
    1.0 0.9 0.6;  % 浅黄色
];

% 绘制旋转平面（每20度一个，共9个）
for i = 0:8
    theta = i * 20 * pi / 180;  % 旋转角度（弧度）
    
    % 旋转矩阵（罗德里格斯公式）
    K = [0, -k(3), k(2); k(3), 0, -k(1); -k(2), k(1), 0];
    R = eye(3) + sin(theta) * K + (1 - cos(theta)) * K^2;
    
    % 旋转平面法向量
    P = [A; B; C];
    P_rotated = R' * P;
    P_rotated = P_rotated / norm(P_rotated);
    A_rot = P_rotated(1);
    B_rot = P_rotated(2);
    C_rot = P_rotated(3);
    
    % 生成平面网格
    [x_plane, y_plane] = meshgrid(linspace(0, 160, 20), linspace(0, 160, 20));
    z_plane = (-A_rot * x_plane - B_rot * y_plane - D) / C_rot;
    
    % 为每个平面分配不同颜色
    color_idx = mod(i, size(colors, 1)) + 1;
    current_color = colors(color_idx, :);
    
    % 绘制平面
    surf(x_plane, y_plane, z_plane, ...
         'FaceAlpha', 0.25, ...
         'FaceColor', current_color, ...
         'EdgeColor', 'none', ...
         'DisplayName', sprintf('平面 %d (%.0f°)', i+1, i*20));
end

% 调整视角和坐标范围
xlim([0 160]);
ylim([0 160]);
zlim([0 160]);

% 添加图例（限制显示的图例数量）
%legend('起始点', '终点', '轴线', 'Location', 'best');

hold off;

% 打印旋转信息
disp('平面旋转完成，共9个平面，每20度一个。');
disp('起始点颜色改为柔和的绿色，轴线使用灰蓝色。');