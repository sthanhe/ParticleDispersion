%% Analyze Stationary Test Measurements
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the paper:
%
%Thanheiser, S.; Haider, M.
%Dispersion Model for Level Control of Bubbling Fluidized Beds with 
%Particle Cross-Flow
%Applied Thermal Energy 2024
%
%All data, along with methodology reports and supplementary documentation, 
%is published in the data repository:
%https://doi.org/10.5281/zenodo.7924693
%
%All required files for this script can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.7948224
%
%
%
%This script analyzes the data of the stationary tests and creates all
%published figures.
%
%
%Requires the file "stat_SumPrep.csv" that summarizes the particle 
%dispersion measurements, which gets created by the script "prepStatic" 
%and stored in the dirStationary folder ("../DataStationary" by default).
%
%Required products:
%   - MATLAB, version 9.14
%   - Simulink, version 10.7
%   - Simulink Real-Time, version 8.2
%   - Stateflow, version 10.8
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @FluBed
%   - @implExp
%   - @Sinter
%   - getBIC.m
%   - loadGeometry.m
%   - postStatic.m
%   - dynamicModel.slx
%   - stat_SumPrep.csv


%% Set data directories
dirStationary='../DataStationary';     %Path to directory where stationary simulation data should be stored
dirFigures='../Figures';               %Path to directory where figures should be stored

%Create directory if it does not exist
if ~isfolder(dirFigures)
    mkdir(dirFigures);
end


%% Load
%Dynamic model, activate fast restart
load_system('dynamicModel');
set_param('dynamicModel',"FastRestart","on");
cleanup=onCleanup(@() set_param('dynamicModel',"FastRestart","off"));

%Data
flow=readtable([dirStationary,filesep,'stat_SumPrep.csv']);
baffleMat=flow{:,compose('baffleCorr%d',1:3)};  %Baffle correction factor matrix

run=1:height(flow);   %Runs to analyze


%% Constant input values
loadGeometry;           %Basic geometry
gates=[false,true];     %Gates (open / close)
direction=true;         %Operating direction, true=forward (left to right)
Y0=ones(1,nACs);        %AC valve rate limiter initial condition


%% Initial state for faster simulations
p0=flow.p0(1);              %Ambient pressure
baffleCorr=baffleMat(1,:);  %Baffle correction factors
Phigate=flow.Phigate(1);    %Weir boundary condition

[bc,Phi,mAC,HAC,mAB]=getBIC(flow(1,:));     %Get other boundary and initial conditions


%Simulate
out=sim('dynamicModel','LoadExternalInput','on','ExternalInput','bc');
xInit=out.xFinal;   %Initial state for other simulations = end state of this simulation


%% Run simulations
for i=run
    p0=flow.p0(i);              %Ambient pressure
    Phigate=flow.Phigate(i);    %Weir boundary condition
    baffleCorr=baffleMat(i,:);  %Baffle correction factors

    bc=getBIC(flow(i,:));   %Get other boundary and initial conditions

    
    %Simulate
    out=sim('dynamicModel','LoadExternalInput','on','ExternalInput','bc',...
                        'LoadInitialState','on','InitialState','xInit');


    %Post processing
    postStatic(out,flow(i,:),x,xChambers,i,dirFigures);
end


%Deactivate fast restart
delete(cleanup);




