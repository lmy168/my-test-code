# PC-RRT*(Planar Constraint RRT*)
The developed planner [PC-RRT*] is an evolution from RRT*.
In most scenarios, the PC-RRT* algorithm exhibits substantial improvements in path length, node efficiency and path smoothness.

# Project Introduction
PC-RRT transforms a 3D path planning problem into a series of 2D path planning problems on multiple planes rotated around an axis determined by the start and goal points.

# Software Requirements
MATLAB R2016a or later.

# File Structure
├── main9.m                   # 主程序入口，障碍物环境1
├── main9zengjiayuanshi.m     # 主程序入口，障碍物环境1，包含未拟合曲线和拟合后的曲线
├── main12.m                  # 主程序入口，障碍物环境2
├── main12zengjiayuanshi.m    # 主程序入口，障碍物环境2，包含未拟合曲线和拟合后的曲线
├── checkpath3.m              # 三维路径碰撞检测（圆形障碍物）
├── feasiblepoint3.m          # 三维路径碰撞检测（方形障碍物）
└── README.md                 # 说明文档

# Quick Start
在安装matlab的windows系统上克隆本项目
在MATLAB R2016a中打开项目文件夹
运行主程序main9或者main12

# Parameter Settings

## 1.基础参数
GoalThreshold = 30;         % 设置目标点阈值
Delta =10;                  % 设置扩展步长
RadiusForNeib = 40;         % rewire的范围，半径r
MaxIterations = 2500;       % 最大迭代次数

## 2.地图参数
searchSize = [250 250 250];  % 搜索范围/地图尺寸

# 致谢


# Maintaince
我们将努力扩展PC-RRT*的适用范围并提高代码的可靠性。

For any technical issues, please contact Wensong Jiang (jwensong@cjlu.edu.cn) or Minyue Li (b24020804011@cjlu.edu.cn).
