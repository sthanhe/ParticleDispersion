%% Baffle Calibration
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
%This script calculates the baffle correction factors to account for their
%influence on particle dispersion. It gets called by either the static or
%dynamic preparation scripts, "prepStatic" or "prepDynamic".
%
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
%   - getBCF.m
%   - getBIC.m
%   - loadGeometry.m
%   - dynamicModel.slx


%% Constant input values
loadGeometry;           %Basic geometry
gates=[false,true];     %Gates (open / close)
direction=true;         %Operating direction, true=forward (left to right)
Y0=ones(1,nACs);        %AC valve rate limiter initial condition


%% Prepare table
hIdx=[1,4,6];   %Bed level indices
names={'PhigateLow','PhigateHigh',...
        'baffleLow','baffleHigh',...
        'hLow','hHigh',...
        'Ylow','Yhigh'};
tab=table('Size',[height(flow),length(names)],...
            'VariableTypes',repmat({'double'},1,length(names)),...
            'VariableNames',names);
clear('names');


%Low and high estimates
tab.PhigateLow=hGate./href.*rho_p.*(1-flow.epsRight)+flow.p0./(FluBed.g.*href);
tab.PhigateHigh=tab.PhigateLow.*1.005;

tab.baffleLow=ones(height(tab),1);
tab.baffleHigh=50*ones(height(tab),1);

ACactive=(flow.AC1~=1)';    %Runs were air cushion is active in the beginning
tab.baffleLow(ACactive)=20;


%Baffle correction factor matrix
baffleMat=ones(height(flow),length(nChambers)-1);


%% Initial state for faster simulations
p0=flow.p0(1);              %Ambient pressure
baffleCorr=baffleMat(1,:);  %Baffle correction factors
Phigate=tab.PhigateLow(1);  %Weir boundary condition

[bc,Phi,mAC,HAC,mAB]=getBIC(flow(1,:));     %Get other boundary and initial conditions


%Simulate
out=sim('dynamicModel','LoadExternalInput','on','ExternalInput','bc');
xInit=out.xFinal;   %Initial state for other simulations = end state of this simulation


%% Phigate (weir boundary condition)
for i=1:height(flow)
    p0=flow.p0(i);          %Ambient pressure
    bc=getBIC(flow(i,:));   %Get other boundary and initial conditions


    %Simulate low Phigate estimate
    Phigate=tab.PhigateLow(i); %#ok<NASGU>

    out=sim('dynamicModel','LoadExternalInput','on','ExternalInput','bc',...
                        'LoadInitialState','on','InitialState','xInit');

    hsim=timeseries2timetable(out.h);
    tab.hLow(i)=hsim.hOverX(end,posBedLevel(end));  %Low h1 estimate


    %Simulate high Phigate estimate
    Phigate=tab.PhigateHigh(i);

    out=sim('dynamicModel','LoadExternalInput','on','ExternalInput','bc',...
                        'LoadInitialState','on','InitialState','xInit');

    hsim=timeseries2timetable(out.h);
    tab.hHigh(i)=hsim.hOverX(end,posBedLevel(end));     %High h1 estimate


    %Get Phigate from estimates through linear interpolation
    flow.Phigate(i)=interp1([tab.hLow(i),tab.hHigh(i)],...
                            [tab.PhigateLow(i),tab.PhigateHigh(i)],...
                            flow.h1(i));
end


%% Baffle 2: inactive AC1
baffleIdx=2; %Baffle index
bafflePos=posBedLevel(3); %#ok<NASGU> %Baffle position index
hname='h4'; %#ok<NASGU> %Name of bed level to be simulated


%Runs to analyze: only where AC1 is inactive in the beginning
run=find(~ACactive);


%Record ACset values, deactivate PIDs
ACsetName=compose('AC%dset',1:nACs);
ACset=flow{run,ACsetName};
flow{run,ACsetName}=ones(length(run),nACs);


%Get baffle correction factors
getBCF;


%Reactivate PIDs
flow{run,ACsetName}=ACset;


%% Baffle 2: active AC1
for i=find(ACactive)
    if isnan(flow.h4(i))
        continue;
    end

    %Boundary conditions
    p0=flow.p0(i);
    baffleCorr=baffleMat(i,:);
    Phigate=flow.Phigate(i);
    bc=getBIC(flow(i,:));


    %Simulate low baffleCorr estimate
    baffleCorr(baffleIdx)=tab.baffleLow(i);

    out=sim('dynamicModel','LoadExternalInput','on','ExternalInput','bc',...
                        'LoadInitialState','on','InitialState','xInit');

    Ysim=timeseries2timetable(out.Ypid);
    tab.Ylow(i)=Ysim.Data(end,1);   %Low AC actuating value estimate


    %Simulate high baffleCorr estimate
    baffleCorr(baffleIdx)=tab.baffleHigh(i);

    out=sim('dynamicModel','LoadExternalInput','on','ExternalInput','bc',...
                        'LoadInitialState','on','InitialState','xInit');

    Ysim=timeseries2timetable(out.Ypid);
    tab.Yhigh(i)=Ysim.Data(end,1);   %High AC actuating value estimate


    %Use bisection method to find baffleCorr
    err=1;
    count=0;
    while err>1e-4 && count<10
        %Simulate new estimate
        baffleCorr(baffleIdx)=interp1([tab.Ylow(i),tab.Yhigh(i)],...
                        [tab.baffleLow(i),tab.baffleHigh(i)],...
                        flow.AC1(i),...
                        'linear','extrap');


        out=sim('dynamicModel','LoadExternalInput','on','ExternalInput','bc',...
                        'LoadInitialState','on','InitialState','xInit');

        Ysim=timeseries2timetable(out.Ypid);
        Ybisect=Ysim.Data(end,1);   %Bisected AC actuating value estimate


        %Get new limits depending on bisected value
        err=abs(Ybisect-flow.AC1(i));
        if Ybisect>flow.AC1(i)
            tab.Ylow(i)=Ybisect;
            tab.baffleLow(i)=baffleCorr(baffleIdx);
        else
            tab.Yhigh(i)=Ybisect;
            tab.baffleHigh(i)=baffleCorr(baffleIdx);
        end


        count=count+1;
    end

    baffleMat(i,:)=baffleCorr;
end


%% Baffle 1
baffleIdx=1;                %Baffle index
bafflePos=posBedLevel(1);   %Baffle position index
hname='h6';                 %Name of bed level to be simulated

tab.baffleLow(ACactive)=1;  %Reset low estimate

run=1:height(flow);     %Runs to analyze


%Get baffle correction factors
getBCF;


%% Add baffle correction matrix to flow table
flow{:,compose('baffleCorr%d',1:length(nChambers)-1)}=baffleMat;




