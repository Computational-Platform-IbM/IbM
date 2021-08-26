% %%%% Script call.m: Initializates the Simulation and calls the numerical
% solver
tic
clc
clear all
% path = [pwd '\shovingQuadTree.jar']; %Windows
path = '../IbM_V2/lib/shovingQuadTree.jar'; %Windows&Cluster
javaaddpath(path,'-end');
% load('R.mat')
R = loadModelXlsx('../IbM_V2/Testing.xlsx');
R.Sxy.pos_xySys = 0;

fprintf('\n> MODEL RUNNING >>>>>\n')
R = integTime(R);
fprintf('\n \n SIMULATION FINISHED >>>> \n')