%% Analyze Mass Diffusivity Measurements
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
%This script analyzes the mass diffusivity calibration data and creates all
%published figures.
%
%
%Requires all mass diffusivity data files ("massDiff_Run...") in the folder
%configured below and all auxiliary classes and functions on the MATLAB 
%path
%
%Required products:
%   - MATLAB, version 9.14
%   - Curve Fitting Toolbox, version 3.9
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @FluBed
%   - @Orifice
%   - @implExp
%   - getConstants.m
%   - getProp.m


%% Set data directories
dirData='DataMassDiff';             %Path to directory containing the data
dirStationary='DataStationary';     %Path to directory where stationary simulation data should be stored
dirFigures='Figures';               %Path to directory where figures should be stored

%Create storage folders if they do not exist
if ~isfolder(dirStationary)
    mkdir(dirStationary);
end

if ~isfolder(dirFigures)
    mkdir(dirFigures);
end


%% Prepare analysis
%Get constants
c=getConstants();


%Retrieve filenames
dirCont=dir(dirData);   %Content of directory containing the data
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'massDiff_'));

%Retain natural order of runs
idx=str2double(extract(files,digitsPattern));
[~,idx]=sort(idx);
files=files(idx);


%Set up table for mean values
chambers=1:6;
names=[{'Run','FG2','Tbed2','wmf2','mDotSand','mDotS','Gamma2','p0','Tleft','Tcenter','Tright','epsLeft','epsCenter','epsRight','AC1','AC2','air1','air2','air3','air4'},...
        strcat(repmat({'h'},1,length(chambers)),arrayfun(@(x) {num2str(x)},chambers)),...
        strcat(repmat({'rho'},1,length(chambers)),arrayfun(@(x) {num2str(x)},chambers)),...
        strcat(repmat({'pi'},1,2),arrayfun(@(x) {num2str(x)},1:2))];
flow=table('Size',[length(files),length(names)],'VariableTypes',repmat({'double'},1,length(names)));
flow.Properties.VariableNames=names;
clear('names');


%% Read individual files and do calculations
for i=1:length(files)
    %Get properties
    tab=readtable([dirData,filesep,files{i}]);
    mDiff=getProp(tab,c,flow.Properties.VariableNames(2:end-2),chambers);   
    
    
    %Calculate bed densities rho=f(mean(eps))
    for j=chambers
        k=num2str(j);

        mDiff{:,['rho',k]}=c.rho_p.*(mDiff{:,['h',k]}+c.hBed)/c.hRef;
        switch j
                case {1,2}
                    mDiff{:,['rho',k]}=mDiff{:,['rho',k]}.*(1-mean(mDiff.epsRight));
                case {3,4}
                    mDiff{:,['rho',k]}=mDiff{:,['rho',k]}.*(1-mean(mDiff.epsCenter));
                case {5,6}
                    mDiff{:,['rho',k]}=mDiff{:,['rho',k]}.*(1-mean(mDiff.epsLeft));
        end
    end
    
    
    %Calculate mass diffusivity
    mDiff.Gamma2=mDiff.mDotS./((mDiff.rho5-mDiff.rho4)./(c.x2-c.xPress));
    
    
    %Remove outliers
    outliers=isoutlier(mDiff.Gamma2,1);
    mDiff{outliers,2:end}=NaN;


    %Record table for future analysis
    writetable(mDiff,[dirStationary,filesep,'stat_Run',num2str(i),'.csv']);
    
    
    %Get means
    flow{i,2:length(mDiff.Properties.VariableNames)}=mean(mDiff{:,2:end},1,'omitnan');
    flow.Run(i)=i;
end


%Add persistent bed levels
for j=chambers
    hloc=['h',num2str(j)];
    flow{:,hloc}=flow{:,hloc}+c.hBed;
end


%Calculate dimensionless variables
flow.pi1=(flow.wmf2.*c.d_p./flow.Gamma2).^-1;
flow.pi2=flow.FG2-1;


%Record table for future analysis
writetable(flow,[dirStationary,filesep,'stat_Sum.csv']);
flow(isnan(flow.Gamma2),:)=[];


%% Fit model to data
%Add boundary condition pi1=pi2=0
pi1=[0;flow.pi1];
pi2=[1e-3;flow.pi2];


%Identify outliers
outliers=[false;...
            (flow.pi2>2 & flow.pi2<2.5 & flow.pi1<53e3) |...
            (flow.pi2>2 & flow.pi2<2.5 & flow.pi1>65e3) |...
            (flow.pi2>2.5 & flow.pi1<68e3)];

%Create figure with outliers
fig=figure(1);
clf(fig);
ax=gca();
box(ax,'on');
hold(ax,'on');

scatter(ax,pi2(~outliers),pi1(~outliers));
scatter(ax,pi2(outliers),pi1(outliers),'x');

hold(ax,'off');

legend({'Valid','Outlier'},'Location','best');

xlabel(ax,'$\pi_2$','Interpreter','latex');
ylabel(ax,'$\pi_1$','Interpreter','latex');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,[dirFigures,filesep,'Outliers.tiff']);


%Get fit and confidence intervals
pi1=pi1(~outliers);
pi2=pi2(~outliers);

[mdl,gof]=fitGamma(pi2,pi1);
ci=confint(mdl);


%Summarize model with fitted coefficients and lower and upper bounds
C=[mdl.a,ci(1,1),ci(2,1)];
a=[mdl.b,ci(1,2),ci(2,2)];

fx=@(x) C(1).*c.hRef/c.hRef0.*x.^a(1);
fxLow=@(x) C(2).*c.hRef/c.hRef0.*x.^a(2);
fxHigh=@(x) C(3).*c.hRef/c.hRef0.*x.^a(3);


%% Create figures
%Figure 5
figidx=5;
fig=figure(figidx);
clf(fig);
ax=gca();
box(ax,'on');
colors=ax.ColorOrder;
hold(ax,'on');


scatter(ax,pi2,pi1);

x=linspace(1e-3,3,100)';
plot(ax,x,fx(x));

plot(ax,x,fxLow(x),'Color',colors(2,:),'LineStyle','--');
plot(ax,x,fxHigh(x),'Color',colors(2,:),'LineStyle','--');

hold off

legend({'Measurements','Fit','95% confidence interval'},'Location','northwest');

xlabel(ax,'\pi_2');
ylabel(ax,'\pi_1');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

name=[dirFigures,filesep,'Figure',num2str(figidx)];
exportgraphics(fig,[name,'.eps']);
exportgraphics(fig,[name,'.tiff']);



%Figure 11
figidx=11;

T=[0;200;400;600]+273.15;
d_p=linspace(80e-6,250e-6,100);
w=4*FluBed.wmf(d_p,c.rho_p,1e5,T);

sz=implExp.size(T,d_p,w);
[Tnorm,d_pNorm,w]=implExp.normalize(sz,T,d_p,w);
gamma=FluBed.Gamma(w,1,d_pNorm,c.rho_p,1e5,Tnorm);
gamma=reshape(gamma,sz);


fig=figure(figidx);
clf(fig);
ax=gca();
box(ax,'on');

plot(ax,d_p.*10^6,gamma);

legend(compose('T=%.0f °C',T-273.15),'Location','northwest');

xlabel(ax,'d_p (\mu m)');
ylabel(ax,'\Gamma (m^2/s)');

ax.XLim=[min(d_p),max(d_p)].*10^6;

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

name=[dirFigures,filesep,'Figure',num2str(figidx)];
exportgraphics(fig,[name,'.eps']);
exportgraphics(fig,[name,'.tiff']);



%Figure 12
figidx=12;

T=[0;200;400;600]+273.15;
d_p=150e-6;
wmf=FluBed.wmf(d_p,c.rho_p,1e5,T);
w=linspace(1.001,6,100).*wmf;

sz=implExp.size(T,d_p,w);
[Tnorm,d_p,wNorm]=implExp.normalize(sz,T,d_p,w);
gamma=FluBed.Gamma(wNorm,1,d_p,c.rho_p,1e5,Tnorm);
gamma=reshape(gamma,sz);


fig=figure(figidx);
clf(fig);
ax=gca();
box(ax,'on');

plot(ax,(w./wmf-1)',gamma');

legend(compose('T=%.0f °C, w_{mf}=%.3f cm/s',T-273.15,wmf.*100),'Location','northwest');

xlabel(ax,'w_e/w_{mf} (-)');
ylabel(ax,'\Gamma (m^2/s)');

ax.XLim=[0,5];

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

name=[dirFigures,filesep,'Figure',num2str(figidx)];
exportgraphics(fig,[name,'.eps']);
exportgraphics(fig,[name,'.tiff']);




%% Fit function
function [fitresult,gof]=fitGamma(pi2,pi1)
%CREATEFIT(PI2,PI1)
%  Create a fit.
%
%  Data for 'mDiffFit' fit:
%      X Input: pi2
%      Y Output: pi1
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 08-May-2023 14:23:51


% Fit: 'mDiffFit'.
[xData, yData] = prepareCurveData( pi2, pi1 );

% Set up fittype and options.
ft = fittype( 'power1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [5435.43821242705 3.67864434815535];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'mDiffFit' );
% h = plot( fitresult, xData, yData );
% legend( h, 'pi1 vs. pi2', 'mDiffFit', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'pi2', 'Interpreter', 'none' );
% ylabel( 'pi1', 'Interpreter', 'none' );
% grid on


end




