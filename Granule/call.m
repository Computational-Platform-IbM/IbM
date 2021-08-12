% %%%% Script call.m: Initializates the Simulation and calls the numerical
% solver
tic
clc
% path = [pwd '\shovingQuadTree.jar']; %Windows
path = [pwd '/shovingQuadTree.jar']; %Windows&Cluster
javaaddpath(path,'-end');
% load('R.mat')
R = loadModelXlsx();
R.Sxy.pos_xySys = 0;

fprintf('\n> MODEL RUNNING >>>>>\n')
R = integTime(R);
fprintf('\n \n SIMULATION FINISHED >>>> \n')