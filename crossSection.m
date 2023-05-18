%% Calculate Relative Free Cross Section
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the paper:
%
%Thanheiser, S.; Haider, M.
%Particle Mass Diffusion Model for Level Control of Bubbling Fluidized Beds
%with Horizontal Particle Flow
%Powder Technology 2023
%
%All data, along with methodology reports and supplementary documentation, 
%is published in the data repository:
%https://doi.org/10.5281/zenodo.7924694
%
%All required files for this script can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.xxxxxxx
%
%
%
%This script calculates the relativ free cross section mentioned in the 
%paper as a possible future parameter in the mass diffusivity model.
%
%
%Requires the file "stat_Sum.csv" that summarizes the mass diffusivity
%measurements, which gets created by the script "calcMassDiff" and stored 
%in the dirStationary folder ("DataStationary" by default).
%
%Required products:
%   - MATLAB, version 9.14
%Necessary files, classes, functions, and scripts:
%   - stat_Sum.csv


%% Set data directory
dirStationary='DataStationary';     %Path to directory where stat_Sum.csv is located


%% Load data
flow=readtable([dirStationary,filesep,'stat_Sum.csv']);


%% Calculations
%Tube geometry
d_tube=20e-3;   %Plain tube diameter
l_tube=500e-3;  %Tube length
n_tube=7;       %Number of tubes in the cross section
l_lead=5e-3;    %Lead distance from the end of the tube to the start of the fins
h_fin=10e-3;    %Fin height
s_fin=2e-3;     %Fin thickness
pitch=9e-3;     %Fin pitch

n_fin=2*floor((l_tube-2*l_lead)/pitch)+1;   %Number of fins
A_fin=h_fin*s_fin*n_fin;                    %Fin cross section area
A_tube=d_tube*l_tube;                       %Plain tube cross section area

A_tubes=n_tube*(A_tube+A_fin);  %Total cross section area of all tubes


h_bed=mean([flow.h4,flow.h5],2);    %Mean bed height of all measurements in the second chamber (only one used for mass diffusivity measurements)

A_0min=l_tube*min(h_bed);       %Minimum bed cross section area without tubes
A_0max=l_tube*max(h_bed);       %Maximum bed cross section area without tubes

CS_freeMin=1-A_tubes/A_0max;    %Minimum relative free cross section area
CS_freeMax=1-A_tubes/A_0min;    %Maximum relative free cross section area




