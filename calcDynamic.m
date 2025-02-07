%% Analyze Dynamic Tests
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the paper:
%
%Thanheiser, S.; Haider, M.
%Dispersion Model for Level Control of Bubbling Fluidized Beds with 
%Particle Cross-Flow
%Chemical Engineering Research and Design 2025
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
%Required products, version 24.1:
%   - MATLAB
%   - Simulink
%   - Requirements Toolbox
%   - Simulink Real-Time
%   - Stateflow
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @FluBed
%   - @implExp
%   - @Sinter
%   - getBIC.m
%   - getConstants.m
%   - mdlPostLoadFx.m
%   - loadGeometry.m
%   - getMdotSstatic.m
%   - postDynamic.m
%   - dynamicModel.slx
%   - dyn_Run...
%   - dyn_Sum.csv
%   - StepResponseFigureInsert.tiff


%% Set data directories
dirStationary='../DataStationary';  %Path to directory where stationary simulation data should be stored
dirData='../DataDynamic';           %Path to directory containing the data
dirFigures='../Figures';            %Path to directory where figures should be stored and where StepResponseFigureInsert.tiff is located


%% Load
%Dynamic model, activate fast restart
mdl='dynamicModel';
sys=load_system(mdl);
mdlPostLoadFx;

set_param(mdl,"FastRestart","on");
cleanup=onCleanup(@() set_param(mdl,"FastRestart","off"));


%Data
flow=readtable([dirData,filesep,'dyn_Sum.csv']);
baffleMat=flow{:,compose('baffleCorr%d',1:length(nChambers)-1)};

run=1:height(flow);     %Runs to analyze


%% Run simulations
const=getConstants();   %Basic constants
for i=run
    %Read table of run
    dyn=readtable([dirData,filesep,'dyn_Run',num2str(i),'.csv']);


    %Initial state
    p0=flow.p0(i);                              %Ambient pressure
    baffleCorr=baffleMat(i,:);                  %Baffle correction factors
    Phigate=flow.Phigate(i);                    %Weir boundary condition

    [bc,Phi,mAC,HAC,mAB]=getBIC(flow(i,:),direction);     %Get other boundary and initial conditions
    
    out=sim(sys,'LoadExternalInput','on','ExternalInput','bc',...
                'LoadInitialState','off');
    xInit=out.xFinal;   %Initial state = end state of stationary simulation


    %Set air cushion setpoint step to measured one
    [~,idx]=find(bc,'Name','hSet');
    bc=setElement(bc,idx,timetable(dyn.Time,dyn{:,compose('AC%dset',1:nACs)}),'hSet');
    

    %Simulate
    out=sim(sys,'LoadExternalInput','on','ExternalInput','bc',...
                'LoadInitialState','on','InitialState','xInit');


    %Post process results
    postDynamic(out,flow(i,:),dyn,const.hBed,posBedLevel,i,dirFigures);
end


%Deactivate fast restart
delete(cleanup);




