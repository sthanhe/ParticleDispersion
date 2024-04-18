%% Analyze Dynamic Tests
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
%This script analyzes the data from the dynamic tests and creates all
%published figures.
%
%
%Requires all dynamic test data files ("dyn_Run...") and the file 
%"dyn_Sum.csv", which get created by the script "prepDynamic" and stored in
%the dirData folder ("../DataDynamic" by default). All auxiliary classes 
%and functions must be on the MATLAB path.
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
%   - getConstants.m
%   - loadGeometry.m
%   - postDynamic.m
%   - dynamicModel.slx
%   - dyn_Run...
%   - dyn_Sum.csv


%% Set data directories
dirStationary='../DataStationary';     %Path to directory where stationary simulation data should be stored
dirData='../DataDynamic';      %Path to directory containing the data
dirFigures='../Figures';       %Path to directory where figures should be stored

%Create storage folder if it does not exist
if ~isfolder(dirFigures)
    mkdir(dirFigures);
end


%% Constant input values
loadGeometry;           %Basic geometry
c=getConstants();       %Basic constants
gates=[false,true];     %Gates (open / close)
direction=true;         %Operating direction, true=forward (left to right)


%% Load
%Dynamic model, activate fast restart
load_system('dynamicModel');
set_param('dynamicModel',"FastRestart","on");
cleanup=onCleanup(@() set_param('dynamicModel',"FastRestart","off"));


%Data
flow=readtable([dirData,filesep,'dyn_Sum.csv']);
baffleMat=flow{:,compose('baffleCorr%d',1:length(nChambers)-1)};

run=1:height(flow);     %Runs to analyze


%% Run simulations
for i=run
    %Read table of run
    dyn=readtable([dirData,filesep,'dyn_Run',num2str(i),'.csv']);


    %Initial state
    Y0=ones(1,nACs); %#ok<PREALL>               %AC valve rate limiter initial condition
    p0=flow.p0(i);                              %Ambient pressure
    baffleCorr=baffleMat(i,:);                  %Baffle correction factors
    Phigate=flow.Phigate(i);                    %Weir boundary condition

    [bc,Phi,mAC,HAC,mAB]=getBIC(flow(i,:));     %Get other boundary and initial conditions
    
    out=sim('dynamicModel','LoadExternalInput','on','ExternalInput','bc');
    xInit=out.xFinal;   %Initial state = end state of stationary simulation


    %Set air cushion setpoint step to measured one
    [~,idx]=find(bc,'Name','hSet');
    bc=setElement(bc,idx,timetable(dyn.Time,dyn{:,compose('AC%dset',1:nACs)}),'hSet');

    %Set PID integrator initial condition to last Ypid value
    Y0=out.Ypid.Data(:,1,end)';
    [~,idx]=find(xInit,'Name','ACpid');
    xInit{idx}.Values.Data=Y0;

    %Simulate
    out=sim('dynamicModel','LoadExternalInput','on','ExternalInput','bc',...
                        'LoadInitialState','on','InitialState','xInit');


    %Post process results
    postDynamic(out,flow(i,:),dyn,c.hBed,posBedLevel,i,dirFigures);
end


%Deactivate fast restart
delete(cleanup);




