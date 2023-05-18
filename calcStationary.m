%% Analyze Stationary Test Measurements
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
%This script analyzes the data of the stationary tests and creates all
%published figures.
%
%
%Requires the file "stat_Sum.csv" that summarizes the mass diffusivity
%measurements, which gets created by the script "calcMassDiff" and stored 
%in the dirStationary folder ("DataStationary" by default).
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
%   - getBIC.m
%   - loadGeometry.m
%   - postStatic.m
%   - stat_Sum.csv


%% Set data directories
dirStationary='DataStationary';     %Path to directory where stationary simulation data should be stored
dirFigures='Figures';               %Path to directory where figures should be stored

%Create directories if they do not exist
if ~isfolder(dirStationary)
    mkdir(dirStationary);
end

if ~isfolder(dirFigures)
    mkdir(dirFigures);
end


%% Load data
flow=readtable([dirStationary,filesep,'stat_Sum.csv']);

run=1:height(flow);   %Runs to analyze

load_system('dynamicModel');


%% Input values
gates=[false,true]; %Gates (open / close)
direction=true;     %Operating direction, true=forward (left to right)


%% Initialize
Y0=ones(1,nACs); 
p0=flow.p0(1);
rhogate=hGate./href.*rho_p.*(1-flow.epsRight(1))+p0./(FluBed.g.*href);     %Sand density at boundary when gate is open
[bc,rho,mAC,HAC,mAB]=getBIC(flow(1,:),1);

out=sim('dynamicModel','LoadExternalInput','on','ExternalInput','bc');
xInit=out.xFinal;


%% Run simulations
for i=run
    p0=flow.p0(i);    %Ambient pressure
    rhogate=hGate./href.*rho_p.*(1-flow.epsRight(i))+p0./(FluBed.g.*href);

    bc=getBIC(flow(i,:),1);
    out=sim('dynamicModel','LoadExternalInput','on','ExternalInput','bc',...
                        'LoadInitialState','on','InitialState','xInit');


    postStatic(out,flow(i,:),x,xChambers,i,dirFigures);
end




