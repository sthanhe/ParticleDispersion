%% Dimensional Analysis
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
%This script sets up the dimensional set shown in the paper.
%
%
%Required products, version 24.1:
%   - MATLAB
%Necessary files, classes, functions, and scripts:
%   - None


%% Dimensional Set
dims={'m','s','kg'};    %Dimensions


%Influencing factors and their dimensions
D=struct('name','D','m',2,'s',-1,'kg',0);
w_e=struct('name','w_e','m',1,'s',-1,'kg',0);
w_p=struct('name','w_p','m',1,'s',-1,'kg',0);
d_p=struct('name','d_p','m',1,'s',0,'kg',0);
rho_p_rho_g=struct('name','rho_p_rho_g','m',-3,'s',0,'kg',1);
rho_g=struct('name','rho_g','m',-3,'s',0,'kg',1);
my_g=struct('name','my_g','m',-1,'s',-1,'kg',1);
g=struct('name','g','m',1,'s',-2,'kg',0);


%Assigning influencing factors to matrices A and B
sA=[d_p,rho_p_rho_g,g];
sB=[D,w_e,w_p,rho_g,my_g];


%Constituting numbers
nA=numel(sA);       %Number of variables in Matrix A
nB=numel(sB);       %Number of variables in Matrix B
nDims=numel(dims);  %Number of dimensions
nP=nA+nB-nDims;     %Number of pi-factors


%Matrices
A=table('Size',[nDims,nA],...
        'VariableTypes',repmat({'double'},1,nA),...
        'VariableNames',{sA.name},...
        'RowNames',dims);

B=table('Size',[nDims,nB],...
        'VariableTypes',repmat({'double'},1,nB),...
        'VariableNames',{sB.name},...
        'RowNames',dims);

C=table('Size',[nP,nA],...
        'VariableTypes',repmat({'double'},1,nA),...
        'VariableNames',{sA.name},...
        'RowNames',compose('pi%d',1:nP));

D=table('Size',[nP,nB],...
        'VariableTypes',repmat({'double'},1,nB),...
        'VariableNames',{sB.name},...
        'RowNames',compose('pi%d',1:nP));


A{:,:}=[sA.m;sA.s;sA.kg];
B{:,:}=[sB.m;sB.s;sB.kg];

D{:,:}=eye(nP);


%Fundamental equation
C{:,:}=-D{:,:}*(A{:,:}^-1*B{:,:})';




