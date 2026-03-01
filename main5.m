k = [1, 1, 1]/sqrt(3); % 空间直线在xyz坐标系下的方向向量

for i = 0:10:360
    theta = i * pi / 180; % 旋转的角度，单位为弧度

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
    fprintf('旋转角度：%d\n', i);
    fprintf('旋转前的平面方程：%.2fx + %.2fy + %.2fz + %.2f = 0\n', A, B, C, D);
    fprintf('旋转后的平面方程：%.2fx + %.2fy + %.2fz + %.2f = 0\n\n', A_rotated, B_rotated, C_rotated, D); 
end