%% Dimensional Analysis
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
%This script checks whether the chosen system of dimensionless parameters
%that make up the particle diffusivity model are linearly independent.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%Necessary files, classes, functions, and scripts:
%   - None


%% Analysis
n=4;        %Number of influencing factors
dims=2;     %Total number of dimensions of the influencing factors


%Matrix of the influencing factors' exponents
A=table('Size',[n,dims],'VariableTypes',repmat({'double'},1,dims),...
        'VariableNames',{'m','s'},...
        'RowNames',{'Gamma','d_p','w_e','w_mf'});

%Mass diffusivity Gamma, m²/s
A.m('Gamma')=2;       
A.s('Gamma')=-1;

%Particle diameter d_p, m
A.m('d_p')=1;

%Excess fluidization velocity w_e=w-w_mf, m/s
A.m('w_e')=1;
A.s('w_e')=-1;

%Minimum fluidization velocity w_mf, m/s
A.m('w_mf')=1;
A.s('w_mf')=-1;


%Matrix of the dimensionless variables' exponents
k=table('Size',[n-dims,n],'VariableTypes',repmat({'double'},1,n),...
        'VariableNames',A.Properties.RowNames);

%pi1=Gamma/(w_mf d_p)
k.Gamma(1)=1;
k.w_mf(1)=-1;
k.d_p(1)=-1;

%pi2=w_e/wmf
k.w_e(2)=1;
k.w_mf(2)=-1;


%Check linear dependence
if all(k{:,:}*A{:,:}==0,'all')
    disp('Variable system is linearly independent');
else
    warning('Variable system is linearly dependent!')
end




