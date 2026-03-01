clc;
%fb = load('C:/Users/15346/Desktop/xyz.txt');
fb = load('C:/Users/15346/Desktop/xyz.txt');
fb = fb(:,1:3);
x = fb(:,1);
y = fb(:,2);
z = fb(:,3);
figure(1)
plot3(x,y,z,'yo'); hold on
xlabel('x');
ylabel('y');
zlabel('z');
% 分别拟合两个二维的曲线，然后统一到一起
p_yx = polyfit(y,x,4);
x_out = polyval(p_yx, y);
p_yz = polyfit(y,z,4);
z_out = polyval(p_yz, y);
plot3(x_out ,y, z_out, 'r*'); hold on;
 
% 得出曲线函数 x_out = f(z_out)   z_out = f(y)
p_zx_out = polyfit(z_out,x_out,4);   
x_out_f = polyval(p_zx_out,z_out);
plot3(x_out_f,y,z_out,'b*'); hold on;