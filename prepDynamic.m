%% Prepare Dynamic Test Analysis
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
%This script prepares the data of the dynamic tests for further analysis
%by the "calcDynamic" script.
%
%
%Requires all dynamic test data files ("dynRaw_Run...") in the folder
%configured below and the file "stat_SumPrep.csv" that summarizes the 
%particle dispersion measurements, which gets created by the script 
%"prepStatic" and stored in the dirStationary folder ("../DataStationary" 
%by default). All auxiliary classes and functions must be on the MATLAB 
%path.
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
%   - @Orifice
%   - baffleCalib.m
%   - getBCF.m
%   - getBIC.m
%   - getConstants.m
%   - getProp.m
%   - mdlPostLoadFx.m
%   - loadGeometry.m
%   - getMdotSstatic.m
%   - dynamicModel.slx
%   - dynRaw_Run...
%   - stat_SumPrep.csv


%% Set data directories
dirStationary='../DataStationary';  %Path to directory where stationary simulation data should be stored
dirData='../DataDynamic';           %Path to directory containing the data
dirFigures='../Figures';            %Path to directory where figures should be stored

%Create storage folder if it does not exist
if ~isfolder(dirFigures)
    mkdir(dirFigures);
end


%% Load dynamic model and activate fast restart
mdl='dynamicModel';
sys=load_system(mdl);
mdlPostLoadFx;

set_param(mdl,"FastRestart","on");
cleanup=onCleanup(@() set_param(mdl,"FastRestart","off"));


%% Prepare analysis
%Get constants
const=getConstants();


%Retrieve filenames
dirCont=dir(dirData);   %Content of directory containing the data
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'dynRaw_'));

%Retain natural order of runs
idx=str2double(extract(files,digitsPattern));
[~,idx]=sort(idx);
files=files(idx);


%Set up table for mean values
flow=readtable([dirStationary,filesep,'stat_SumPrep.csv']);
flow{:,2:end}=0;
flow(length(files)+1:end,:)=[];
names=flow.Properties.VariableNames(2:end);
chambers=1:6;


%% Read individual files and do calculations
%Variables whose stationary value in the beginning is needed
idxStat=[compose('AC%d',1:2),compose('AC%dset',1:2),...
        compose('h%d',1:6),compose('Phi%d',1:6)];

for i=1:length(files)
    %Read file
    tab=readtable([dirData,filesep,files{i}]);

    %Set up table
    % names=[{'Time'},flow.Properties.VariableNames(2:end)];
    % dyn=table('Size',[height(tab),length(names)],'VariableTypes',[{'duration'},repmat({'double'},1,length(names)-1)]);
    % dyn.Properties.VariableNames=names;
    % clear('names');

    %Get properties
    % dyn(:,1:end)=getProp(tab,c,dyn.Properties.VariableNames(2:end),chambers);
    dyn=getProp(tab,const,names,chambers);


    %Record table for future analysis
    writetable(dyn,[dirData,filesep,'dyn_Run',num2str(i),'.csv']);
    
    
    %Get means
    flow{i,2:end}=mean(dyn{:,2:end},1,'omitnan');
    flow.Run(i)=i;

    stat=dyn.AC1set==dyn.AC1set(1);     %Identify stationary beginning
    flow{i,idxStat}=mean(dyn{stat,idxStat},1,'omitnan');
end


%Add persistent bed levels
hNames=compose('h%d',chambers);
flow{:,hNames}=flow{:,hNames}+const.hBed;


%% Do baffle calibration
baffleCalib;


%Record table for future analysis
writetable(flow,[dirData,filesep,'dyn_Sum.csv']);




