%% Dynamic Model Post Load Function
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the paper:
%
%Thanheiser, S.; Haider, M.
%Dispersion Model for Level Control of Bubbling Fluidized Beds with 
%Particle Cross-Flow
%Chemical Engineering Science 2024
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
%This script creates the workspace variables and an initial set of 
%boundary and initial conditions necessary to run the dynamic numerical 
%model "dynamicModel.slx". 
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @FluBed
%   - @implExp
%   - @Sinter
%   - getBIC.m
%   - loadGeometry.m
%   - getMdotSstatic.m


%% Basic geometry
loadGeometry;


%% Default values
direction=true;         %Flow direction
baffleCorr=ones(1,3);   %Baffle correction factors


%% Initial set of boundary and inital conditions
%Set up table
names={'p0','Tleft','Tcenter','Tright','epsLeft','epsCenter','epsRight','mDotS','air1','air2','air3','air4','AC1set','AC2set'};
flow=table('Size',[1,length(names)],'VariableTypes',[repmat({'double'},1,length(names))]);
flow.Properties.VariableNames=names;
clear('names');


%Take values directly from a measurement point
flow.p0=101322.321749582;
flow.Tleft=321.622691993331;
flow.Tcenter=319.331017705775;
flow.Tright=330.764596384106;
flow.epsLeft=0.467551270924290;
flow.epsCenter=0.470101396127068;
flow.epsRight=0.472651521329844;
flow.mDotS=4;
flow.air1=0.0141454434627567;
flow.air2=0.0425907793606724;
flow.air3=0.0370863059211219;
flow.air4=0.0121968546724413;
flow.AC1set=1;
flow.AC2set=1;


%Get boundary and initial conditions
[bc,Phi,mAC,HAC,mAB]=getBIC(flow(1,:),direction);


%Other boundary and initial conditions
p0=flow.p0(1);  %Ambient pressure
Phigate=hGate./href.*rho_p.*(1-flow.epsRight(1))+p0./(FluBed.g.*href);  %Weir boundary condition
Y0=ones(1,nACs);    %PID I-value


%% Simulation time
minSimTime=240;      %Minimum simulation time (seconds)
maxSimTime=600;    %Maximum simulation time (seconds)
statCond=1e-4;      %Condition for stationary status as a value of PhiDot


%% Static simulation parameters
mDotSstatic=getMdotSstatic(flow.mDotS,direction);
isStatic=false;

YMan=[1,1];
SetMan=false;


%% Particle dispersion coefficients
c=19030.1124724027;
eps2=1.08175825313772;
eps3=-3.62266055421733;
epsAr=0.109738530498836;




