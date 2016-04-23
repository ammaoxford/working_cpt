clear all; close all; clc

addpath('mfiles/')
addpath('data/mat_Y')
datanames_Y;
load(strcat('mat_Y/',dat_Y{18}))

% t  = t(1:3000);
% bp = bp(1:3000);



%% COMPUTE BRS
[brs,vv,brsdata] = cpt_brs(t,bp);
brs

