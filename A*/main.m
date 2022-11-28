% Path planning
% Thanks to HKUST ELEC 5660 
close all; clear all; clc;

% 2D grid map
xStart = 1.0;
yStart = 1.0;
xTarget = 20.0;
yTarget = 10.0;
MAX_X = 20;
MAX_Y = 10;
map = obstacle_map(xStart, yStart, xTarget, yTarget, MAX_X, MAX_Y);

% Waypoint Generator Using the A* 
path = A_star_search(map, MAX_X,MAX_Y);

% visualize the 2D grid map
visualize_map(map, path);
