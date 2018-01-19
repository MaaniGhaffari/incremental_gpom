clc; clear; close all

% Initialize GPML
addpath('gpml-matlab-v4.0')
startup

load('gpom_test')

output = IGPOM(robotPose, laserScan, Parameters, gpom_store)