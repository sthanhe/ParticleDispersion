%% Post Processing of Stationary Simulations
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
%All required files for this function can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.7948224
%
%
%
%This function anaylzes the results of the stationary simulations conducted
%with the script "calcStationary" and creates all published figures.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary files, classes, functions, and scripts:
%   - None


function postStatic(out,flow,x,xChambers,figidx,dirFigures)
    %% Common parameters
    xmeas=[181.5e-3,37e-3,1023e-3,37e-3,763e-3,37e-3];
    xmeas=cumsum(xmeas);    %x-coordinates of bed level measurements


    %Figure title:
    titleText=['Test ',num2str(figidx),...
                ', $w_{e}$/$w_{mf}$=',num2str(round(flow.FG2-1,1)),...
                ', $\dot{m}_S$=',num2str(flow.mDotS),' kg/m$^{2}$s',...
                ', T=',num2str(round(flow.Tbed2-273.15)),'$^{\circ}$C'];


    %% Bed level
    hmeas=flow{:,compose('h%d',6:-1:1)};    %Measured bed levels
    hsim=timeseries2timetable(out.h);       %Simulated bed levels


    %Set up figure
    fig912=figure(figidx);
    clf(fig912);
    ax=gca();
    box(ax,'on');
    hold(ax,'on');
    colors=ax.ColorOrder;

    plot(ax,x,hsim{end,:},'Color',colors(1,:));
    plot(ax,xmeas,hmeas,'Color',colors(2,:));


    %Baffle positions
    xline(cumsum(xChambers(1:end-1)));
    hold(ax,'off');

    
    %Configure and save figure
    legend(ax,{'Simulated','Measured'});

    title(ax,titleText,'Interpreter','latex');
    xlabel(ax,'Distance from Inlet (m)');
    ylabel(ax,'Bed Level (m)');

    fig912.Units='centimeters';
    fig912.Position=[0.02,12.18,17,8.5];
    ax.XLim=[0,max(x)];

    exportgraphics(fig912,[dirFigures,filesep,'statLevels',num2str(figidx),'.tiff']);


    %% Fluidization
    FG=timeseries2timetable(out.FG);    %Degree of fluidizations simulated
    D=timeseries2timetable(out.D);      %Mass diffusivity simulated


    %Set up figure
    fig913=figure(figidx+100);
    clf(fig913);
    ax=gca();
    box(ax,'on');


    %Excess fluidization
    plot(ax,x,FG{end,:}-1)
    ylabel(ax,'w_e/w_{mf} (-)');


    %Particle dispersion
    hold(ax,'on');
    yyaxis(ax,'right');
    plot(ax,x,squeeze(D{end,:}))
    ylabel(ax,'D (m^2/s)');
    

    %Baffle positions
    xline(cumsum(xChambers(1:end-1)));
    hold(ax,'off');
    

    %Configure and save figure
    title(ax,titleText,'Interpreter','latex');
    xlabel(ax,'Distance from Inlet (m)');

    fig913.Units='centimeters';
    fig913.Position=[0.02,0.83,17,8.5];
    ax.XLim=[0,max(x)];

    exportgraphics(fig913,[dirFigures,filesep,'statFluidization',num2str(figidx),'.tiff']);
end




