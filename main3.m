% 平面方程的系数
A = 1;
B = -1;
C = 0;
D = 0;

% 点的坐标
x0 = 1;
y0 = -2;
z0 = 5;

% 判断点是否在平面上
onPlane = isPointOnPlane(A, B, C, D, x0, y0, z0);

% 输出结果
if onPlane
    disp('点 P 在平面上。');
else
    disp('点 P 不在平面上。');
end
