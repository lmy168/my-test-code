# PC-RRT*(Planar Constraint RRT*)
The proposed planner, PC-RRT*, is an evolution of the RRT algorithm.
In most scenarios, the PC-RRT* algorithm exhibits substantial improvements in path length, node efficiency, and path smoothness.

# Project Introduction
PC-RRT transforms a 3D path planning problem into a series of planar path planning problems on multiple planes rotated around an axis determined by the start and goal points.

# Software Requirements
MATLAB R2016a or later.

# File Structure
|Filename|Description|
|-|-|
|main9.m|Main program entry for obstacle environment 1|
|main9_1.m|Main entry for environment 1 (unfitted and fitted curves in separate figures)|
|main9_2.m|Main entry for environment 1 (unfitted and fitted curves in the same figure)|
|main12.m|Main program entry for obstacle environment 2|
|main12_1.m|Main entry for environment 2 (unfitted and fitted curves in separate figures)|
|main12_2.m|Main entry for environment 2 (unfitted and fitted curves in the same figure)|
|checkpath3.m|3D path collision detection (spherical obstacles)|
|feasiblepoint3.m|3D path collision detection (cuboid obstacles)|
|README.md|Documentation|

# Quick Start
Clone this repository on a Windows system that has MATLAB installed.
Open this project folder with MATLAB R2016a.
Execute the main script main9 or main12.

# Parameter Settings

## 1. Basic Parameters
GoalThreshold = 30;          % Set goal threshold
Delta = 10;                   % Set expansion step size
RadiusForNeib = 40;          % Rewiring range, radius r
MaxIterations = 2500;        % Maximum number of iterations

## 2. Map Parameters
searchSize = [250 250 250];  % 3D map dimensions [x y z]

# License
This project is licensed under the MIT License - see the [MIT License](LICENSE) file for details.

# Maintenance
We aim to generalize PC-RRT to broader applications and further improve its code robustness in future work.

For any technical issues, please contact Wensong Jiang (jwensong@cjlu.edu.cn) or Minyue Li (b24020804011@cjlu.edu.cn).
