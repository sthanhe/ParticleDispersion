%% Operating Point Test
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
%This script can be used to test whether an operating point created by the
%"calcRobust" script is actually in stationary conditions. It requires the
%operating points stored in the "oppoints.mat" file created by the script
%"calcRobust".
%
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
%   - mdlPostLoadFx.m
%   - loadGeometry.m
%   - getMdotSstatic.m
%   - dynamicModel.slx
%   - oppoints.mat


%% Set point to check
opIdx=23;       %Index in the operating point grid
paramIdx=5;     %Index in the parameter grid


%% Load
%Dynamic model
mdl='dynamicModel';
sys=load_system(mdl);
mdlPostLoadFx;


%Operating points and parameter grids
load('oppoints.mat');


%% Check operating point
%Boundary conditions
poro=0.47;  %Bed porosity (constant)
p0=101325;  %Ambient pressure (constant)

Phigate=hGate./href.*rho_p.*(1-poro)+p0./(FluBed.g.*href);  %Weir boundary condition


%Get dispersion model parameters
c=cGrid(paramIdx);
eps2=eps2Grid(paramIdx);
eps3=eps3Grid(paramIdx);
epsAr=epsArGrid(paramIdx);


%Get operating point boundary and initial conditions
point=op{opIdx}(paramIdx);
uInit=getinputstruct(point);
xInit=getstatestruct(point);


%Simulate operating point
minSimTime=100;     %Ensure short simulation
out=sim(sys,'LoadExternalInput','on','ExternalInput','uInit',...
            'LoadInitialState','on','InitialState','xInit');




