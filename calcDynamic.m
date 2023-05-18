%% Analyze Dynamic Tests
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
%This script analyzes the data from the dynamic tests and creates all
%published figures.
%
%
%Requires all dynamic test data files ("dynRaw_Run...") in the folder
%configured below and all auxiliary classes and functions on the MATLAB 
%path
%
%Required products:
%   - MATLAB, version 9.14
%   - Simulink, version 10.7
%   - Simulink Real-Time, version 8.2
%   - Stateflow, version 10.8
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @FluBed
%   - @Orifice
%   - @implExp
%   - getBIC.m
%   - getConstants.m
%   - getProp.m
%   - loadGeometry.m
%   - postDynamic.m


%% Set data directories
dirData='DataDynamic';      %Path to directory containing the data
dirFigures='Figures';       %Path to directory where figures should be stored

%Create storage folder if it does not exist
if ~isfolder(dirFigures)
    mkdir(dirFigures);
end


%% Prepare analysis
%Get constants
c=getConstants();


%Retrieve filenames
dirCont=dir(dirData);   %Content of directory containing the data
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'dynRaw_'));

%Retain natural order of runs
idx=str2double(extract(files,digitsPattern));
[~,idx]=sort(idx);
files=files(idx);


%Set up table for mean values
chambers=1:6;
names=[{'Run','FG2','Tbed2','wmf2','mDotSand','mDotS','p0','Tleft','Tcenter','Tright','epsLeft','epsCenter','epsRight','AC1','AC2','air1','air2','air3','air4'},...
        strcat(repmat({'h'},1,length(chambers)),arrayfun(@(x) {num2str(x)},chambers))];
flow=table('Size',[length(files),length(names)],'VariableTypes',repmat({'double'},1,length(names)));
flow.Properties.VariableNames=names;
clear('names');


%% Read individual files and do calculations
for i=1:length(files)
    %Read file
    tab=readtable([dirData,filesep,files{i}]);

    %Set up table
    names=[{'Time'},flow.Properties.VariableNames(2:end),{'AC1set','AC2set'}];
    dyn=table('Size',[height(tab),length(names)],'VariableTypes',[{'duration'},repmat({'double'},1,length(names)-1)]);
    dyn.Properties.VariableNames=names;
    clear('names');

    %Get properties
    dyn(:,1:end-2)=getProp(tab,c,dyn.Properties.VariableNames(2:end-2),chambers);


    %Write air cushion setpoints
    dyn.AC1set=tab.ACset1+c.hBed;
    dyn.AC2set=tab.ACset2+c.hBed;


    %Record table for future analysis
    writetable(dyn,[dirData,filesep,'dyn_Run',num2str(i),'.csv']);
    
    
    %Get means
    flow{i,2:length(dyn.Properties.VariableNames)-2}=mean(dyn{:,2:end-2},1,'omitnan');
    flow.Run(i)=i;
end


%Add persistent bed levels
for j=chambers
    hloc=['h',num2str(j)];
    flow{:,hloc}=flow{:,hloc}+c.hBed;
end


%Record table for future analysis
writetable(flow,[dirData,filesep,'dyn_Sum.csv']);


%% Load data for dynamic analysis
run=1:height(flow);   %Runs to analyze

load_system('dynamicModel');


%% Input values
gates=[false,true]; %Gates (open / close)
direction=true;     %Operating direction, true=forward (left to right)

%Baffle correction factors
bcf=ones(height(flow),2);
bcf(1,1)=23;
bcf(1,2)=74;
bcf(2,1)=17;
bcf(2,2)=64;
bcf(3,1)=12;
bcf(3,2)=68;
bcf(4,1)=15;
bcf(4,2)=64;
bcf(5,1)=22;
bcf(5,2)=73;
bcf(6,1)=60;
bcf(6,2)=129;
bcf(7,1)=60;
bcf(7,2)=128;
bcf(8,1)=50;
bcf(8,2)=125;
bcf(9,1)=69;
bcf(9,2)=138;
bcf(10,1)=69;
bcf(10,2)=128;
bcf(11,1)=65;
bcf(11,2)=130;


%% Run simulations
for i=run
    %Read table of run, set baffle correction factor
    dyn=readtable([dirData,filesep,'dyn_Run',num2str(i),'.csv']);
    baffleCorr=bcf(i,:);


    %Initialize
    Y0=ones(1,nACs);                                                        %#ok<PREALL> %Rate limiter initial condition
    p0=flow.p0(i);                                                          %Ambient pressure
    rhogate=hGate./href.*rho_p.*(1-flow.epsRight(i))+p0./(FluBed.g.*href);  %Sand density at boundary when gate is open
    [bc,rho,mAC,HAC,mAB]=getBIC(flow(i,:),dyn.AC1set(1));                   %Boundary and initial conditions
    
    out=sim('dynamicModel','LoadExternalInput','on','ExternalInput','bc');
    xInit=out.xFinal;   %Initial conditions for actual simulation


    %Set air cushion setpoint step to measured one
    [~,idx]=find(bc,'Name','hSet');
    bc=setElement(bc,idx,timetable(dyn.Time,dyn.AC1set),'hSet');

    %Set PID integrator initial condition to last Ypid value
    Y0=out.Ypid.Data(1,:,end);          %Rate limiter initial condition
    [~,idx]=find(xInit,'Name','ACpid');
    xInit{idx}.Values.Data=Y0;

    %Simulate
    out=sim('dynamicModel','LoadExternalInput','on','ExternalInput','bc',...
                        'LoadInitialState','on','InitialState','xInit');


    %Post process results
    postDynamic(out,flow(i,:),dyn,c.hBed,posBedLevel,i,dirFigures);
end




