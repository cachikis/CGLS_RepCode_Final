%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: main_cgls.m
% Author: Craig A. Chikis
% Date: 01/06/2022
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -global
clear all
close all

% Fill this in per your local file setup
cd("../../")


addpath('Code/Main/')
addpath('Code/Model/')
addpath('Code/Analysis/')
addpath('Code/Parameterizations/')
addpath('Code/Data/')

tic
run_paper(); 
toc

