%% Analyze Controller Robustness
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
%This script linearizes the dynamic model at different operating points and
%analyzes controller robustness in the frequency domain.
%
%
%Required products, version 24.1:
%   - MATLAB
%   - Robust Control Toolbox
%   - Simulink Control Design
%   - Deep Learning HDL Toolbox
%   - HDL Verifier
%   - Fixed-Point Designer
%   - MATLAB Coder
%   - Signal Processing Toolbox
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
%   - dynamicModel.slx


%% Set data directories
dirTemp='Temp';                        %Path to temporary file directory
dirFigures='../Figures';               %Path to directory where figures should be stored

%Create directories if they do not exist
if ~isfolder(dirTemp)
    mkdir(dirTemp);
end

if ~isfolder(dirFigures)
    mkdir(dirFigures);
end


%% Inputs
mdl='dynamicModel';     %Model name


Tbed=[50,150]+273.15;   %Bed temperature range
mDotS=[4,6];            %Specific particle flow range

FG=[2.5,4];     %Degree of fluidization range
FGhigh=3;       %Degree of fluidization above which YManHigh is used

YManLow=[0.3,0.9];      %AC actuating value range for low FGs
YManHigh=[0.5,0.9];     %AC actuating value range for high FGs

dirGrid=[true,false];   %Particle flow direction range


poro=0.47;      %Bed porosity (constant)
p0act=101325;   %Ambient pressure (constant)


%Standard errors of particle dispersion coefficients and exponents
c_SE=3945.22646777468;
eps2_SE=0.0475987478498504;
eps3_SE=0.485252232546785;
epsAr_SE=0.0371440626730409;


%% Storage management
%Set up temp directory to collect only relevant temp files
clear('tempdir');

tempDirNew=[cd,filesep,dirTemp];
if ispc
    setenv('TMP',tempDirNew)
else
    setenv('TMPDIR',tempDirNew)
end
clear('tempDirNew');

tempdir; %refresh


%Make sure that clearing the temp folder permanently deletes files and
%throws no warnings regarding permissions
recycle('off');
warning('off','MATLAB:DELETE:Permission');


%Deactivate Simulink auto archiving to save disk space
Simulink.sdi.setAutoArchiveMode(false);


%Reverse changes on cleanup (temp directory gets restored after Matlab restarts)
cleanup=onCleanup(@() warning('on','MATLAB:DELETE:Permission'));
cleanup=[cleanup,onCleanup(@(x) Simulink.sdi.setAutoArchiveMode(true))]; %#ok<NASGU>


%% Set up boundary conditions
%Load dynamic model
sys=load_system(mdl);
mdlPostLoadFx;


%Boundary conditions
flow.p0=p0act;
flow.epsLeft=poro;
flow.epsCenter=poro;
flow.epsRight=poro;
flow.AC1set=1;
flow.AC2set=1;

p0=p0act;
Phigate=hGate./href.*rho_p.*(1-poro)+p0./(FluBed.g.*href);  %Weir boundary condition


%Set up stationary calculations
SetMan=true;    %#ok<NASGU>
isStatic=true;  %#ok<NASGU>
minSimTime=0;   %Minimize calculation time


%% Set initial conditions
out=sim(sys,'LoadExternalInput','on','ExternalInput','bc',...
            'LoadInitialState','off');

Phi=get(out.xFinal,'Phi').Values.Data;
mAC=get(out.xFinal,'mAC').Values.Data;
HAC=get(out.xFinal,'HAC').Values.Data;
mAB=get(out.xFinal,'mAB').Values.Data;


%% Operating point and parameter grids
%Operating points
[YMan1,YMan2,Tbed,mDotS,FG,dirGrid]=...
    ndgrid(YManLow,YManLow,Tbed,mDotS,FG,dirGrid);

isFGhigh=FG>FGhigh;
for i=1:length(YManLow)
    YMan1(isFGhigh & YMan1==YManLow(i))=YManHigh(i);    
    YMan2(isFGhigh & YMan2==YManLow(i))=YManHigh(i);
end


%Parameters
cGrid=[c-c_SE,c+c_SE];
eps2Grid=[eps2-eps2_SE,eps2+eps2_SE];
eps3Grid=[eps3-eps3_SE,eps3+eps3_SE];
epsArGrid=[epsAr-epsAr_SE,epsAr+epsAr_SE];

[cGrid,eps2Grid,eps3Grid,epsArGrid]=...
    ndgrid(cGrid,eps2Grid,eps3Grid,epsArGrid);


param=repmat(struct('Name','c','Value',cGrid),1,4);
param(2).Name='eps2';
param(3).Name='eps3';
param(4).Name='epsAr';

param(2).Value=eps2Grid;
param(3).Value=eps3Grid;
param(4).Value=epsArGrid;


%% Calculate operating point snapshots
const=getConstants();
pBed=p0+FluBed.deltaP(const.hBed./2,poro,rho_p);    %Approximate mean bed pressure for all operating points


op=cell(1,numel(Tbed));
for i=1:length(op)
    %Set boundary conditions
    flow.Tleft=Tbed(i);
    flow.Tcenter=Tbed(i);
    flow.Tright=Tbed(i);
    flow.mDotS=mDotS(i);
    
    
    rho_g=DryAir.rho(pBed,Tbed(i));             %Fluidization gas density
    wmf=FluBed.wmf(d_p,rho_p,pBed,Tbed(i));     %Minimum fluidization velocity
    w_e=(FG(i)-1).*wmf;                         %Excess fluidization velocity
    
    flow.air1=(w_e+wmf).*(rho_g.*const.Afloor1);
    flow.air2=(w_e+wmf).*(rho_g.*const.Afloor2);
    flow.air3=(w_e+wmf).*(rho_g.*const.Afloor3);
    flow.air4=flow.air1;
    
    
    bc=getBIC(flow,dirGrid(i));


    %Stationary mDotS and YMan values
    mDotSstatic=getMdotSstatic(flow.mDotS,dirGrid(i));
    YMan=[YMan1(i),YMan2(i)];


    %Get operating points
    points=findop(mdl,1e4,param);
    stateNames={points(1).States.StateName};
    inputNames={points(1).Inputs.Block};


    %Set PID integrator states to YMan values
    idx=find(strcmp(stateNames,'ACpid'));
    for j=1:numel(points)
        points(j).States(idx).x=YMan';
    end


    %Retrieve measured bed levels at stationary states
    idx=find(strcmp(stateNames,'hSave'));
    hout=cell2mat(...
            arrayfun(@(i) ...
                points(i).States(idx).x',...
            1:numel(points),'UniformOutput',false)');


    %Set PID setpoints (inports) to measured bed levels
    idx=find(contains(inputNames,'/hSet'));
    for j=1:numel(points)
        points(j).Inputs(idx).u=hout(j,:)';
    end


    %Set mDotS unit delay to mDotSstatic
    idx=find(strcmp(stateNames,'mDotS'));
    for j=1:numel(points)
        points(j).States(idx).x=mDotSstatic';
    end


    %Save operating points in cell array
    op{i}=points;


    %Clear temp folder
    delete([dirTemp,filesep,'*']);
end


%Save for future analysis
save('oppoints',"op","cGrid","eps2Grid","eps3Grid","epsArGrid");


%% Expand operating point and parameter grids
cGrid=repmat(cGrid,length(op),1);
eps2Grid=repmat(eps2Grid,length(op),1);
eps3Grid=repmat(eps3Grid,length(op),1);
epsArGrid=repmat(epsArGrid,length(op),1);

param(strcmp({param.Name},'c')).Value=cGrid;
param(strcmp({param.Name},'eps2')).Value=eps2Grid;
param(strcmp({param.Name},'eps3')).Value=eps3Grid;
param(strcmp({param.Name},'epsAr')).Value=epsArGrid;


op=vertcat(op{:});

sz=[1,size(cGrid,2:ndims(cGrid))];
YMan1=getGrid(YMan1,sz);
YMan2=getGrid(YMan2,sz);
Tbed=getGrid(Tbed,sz);
mDotS=getGrid(mDotS,sz);
FG=getGrid(FG,sz);
dirGrid=getGrid(dirGrid,sz);


%% Check if all operating points are fluidized
names={op(1).States.StateName};
idx=find(strcmp(names,'FGsave'));
isFluidized=all(...
                arrayfun(@(i) ...
                    all(op(i).States(idx).x>0),...
                1:numel(op)));


%% Linearizations
%Deactivate stationary calculations
SetMan=false;
isStatic=false;


%Controller
ios(1)=linio([mdl,'/error'],1,'openinput');
ios(2)=linio([mdl,'/Controller'],1,'output');

contr=linearize(mdl,ios,op,param);
contr=reshape(contr,numel(op),1);


%Plant
ios(1)=linio([mdl,'/Controller'],1,'openinput');
ios(2)=linio([mdl,'/Plant'],4,'output');

plant=linearize(mdl,ios,op,param);
plant=reshape(plant,numel(op),1);


%Open loop response
openLoop=plant*contr;
Ts=openLoop(1,1,1).Ts;  %Sample time


%Others
% closedLoop=feedback(openLoop,eye(size(openLoop,1)));  %Closed loop
% sensOut=feedback(eye(size(openLoop,1)),openLoop);     %Output sensitivity at bed level measurements


%% Open loop bode plot
%Frequency response at Nyquist frequency
clear('i');             %clear complex variable
w=pi/Ts;                %rad/s at Nyquist frequency         
z=exp(i*w*Ts);          %Laplace-variable


mag=NaN(numel(op),size(openLoop,1));
for j=1:size(openLoop,1)
    mag(:,j)=squeeze(evalfr(openLoop(j,j),z));
end

mag=mag2db(abs(mag));


%Critical magnitudes
idxMin=mag(:,1)<-36 | mag(:,2)<-40;
idxMax=mag(:,1)>-7 | mag(:,2)>-10;


%Set up figure
fig=figure(912);
clf(fig);
ax=gca();


opts=bodeoptions;
opts.FreqUnits='Hz';
opts.PhaseVisible='off';
opts.Title.String='';
opts.Grid='on';
opts.XLimMode='manual';
opts.XLim={[1e-6,1./(2.*Ts)]};

bodeplot(ax,openLoop(:,:,idxMin),openLoop(:,:,idxMax),opts);


legend(ax,{'Min','Max'},'Location','best');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,[dirFigures,filesep,'bodeplot.tiff'],...
    'Resolution',600);


%% Disk margins
%Multiloop margin at plant output
[~,MM]=diskmargin(openLoop);
MM=reshape(MM,[],1);


%Lowest disk margin: YMan1=low, Tbed=low, mDotS=high, dir=reverse, eps3=low
%Highest disk margin: YMan1=high, YMan2=high, Tbed=high, FG=low
DM=[MM.DiskMargin];
idxMin=DM<0.8;
idxMax=DM>1.9;


%Set up figure
fig=figure(913);
clf(fig);
ax=gca();


opts=diskmarginoptions;
opts.FreqUnits='Hz';
opts.Title.String='';
opts.XLimMode='manual';
opts.XLim=[1e-4,1./(2.*Ts)];

diskmarginplot(ax,openLoop(:,:,idxMax),openLoop(:,:,idxMin),opts);


legend(ax,{'Max','Min'},'Location','best');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,[dirFigures,filesep,'diskmargins.tiff'],...
    'Resolution',600);


%% Cleanup
clear('cleanup');


%% Auxiliary function
function x=getGrid(x,sz)
    x=cell2mat(...
        arrayfun(@(i) ...
            repmat(x(i),1,2),...
        1:numel(x),'UniformOutput',false))';

    x=repmat(x,sz);
end




