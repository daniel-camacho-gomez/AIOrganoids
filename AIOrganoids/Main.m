%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Code for the morphogenesis of solid organoids
% Agent-based model (lattice-free center-based) + Neural Network
% By Daniel Camacho-Gomez,
% Unversity of zaragoza, Spain
% E-mail: dcamacho@unizar.es, danielcamachogomez@hotmail.com
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
clear all; close all

%Load the weights of the trained Neural Network
addpath(genpath('Data'))
load('NN_weights.mat') 

%Target Data to simulate 
%...target times
t_t    = [3 5 7]; %days
%...target number of cells for the different target times
n_ct   = [6 20 60]; %number of cells


%Call the agent-based model 
 [n_c_evo] = aiOrganoids(t_t, n_ct, h_1, h_2, o_1, o_2);
                             







