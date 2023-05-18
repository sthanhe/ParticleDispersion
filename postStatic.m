%% Post Processing of Stationary Simulations
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
%All required files for this function can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.xxxxxxx
%
%
%
%This function anaylzes the results of the stationary simulations conducted
%with the script "calcStationary" and creates all published figures.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%Necessary files, classes, functions, and scripts:
%   - None


function postStatic(out,flow,x,xChambers,figidx,dirFigures)
    %% Common parameters
    xmeas=[181.5e-3,37e-3,1023e-3,37e-3,763e-3,37e-3];
    xmeas=cumsum(xmeas);    %x-coordinates of bed level measurements

    baffle=63;  %Index of baffle position

    %Figure title:
    titleText=['Test ',num2str(figidx),...
                ', $w_{e}$/$w_{mf}$=',num2str(round(flow.FG2-1,1)),...
                ', $\dot{m}_S$=',num2str(flow.mDotS),' kg/m$^{2}$s',...
                ', T=',num2str(round(flow.Tbed2-273.15)),'$^{\circ}$C'];


    %% Bed level
    hmeas=[flow.h6,flow.h5,flow.h4,flow.h3,flow.h2,flow.h1];    %Measured bed levels
    hsim=timeseries2timetable(out.h);                           %Simulated bed levels


    %Set up figure
    fig6=figure(figidx);
    clf(fig6);
    ax=gca();
    box(ax,'on');
    hold(ax,'on');
    colors=ax.ColorOrder;

    %Without correction
    plot(ax,x,hsim{end,:},'Color',colors(1,:));
    plot(ax,xmeas,hmeas,'Color',colors(2,:));
    

    %With correction
    hcorr=hmeas;
    hcorr(1:3)=hcorr(1:3)-(hcorr(3)-hsim{end,1}(baffle));
    hcorr(4:end)=hcorr(4:end)-(hcorr(4)-hsim{end,1}(baffle+1));

    plot(ax,xmeas,hcorr,'Color',colors(2,:),'LineStyle','--');


    %Baffle positions
    xline(cumsum(xChambers(1:end-1)));
    hold(ax,'off');

    
    %Configure and save figure
    legend(ax,{'Simulated','Measured','Corrected'});

    t6=title(ax,titleText,'Interpreter','latex');
    xlabel(ax,'Distance from Inlet (m)');
    ylabel(ax,'Bed Level (m)');

    fig6.Units='centimeters';
    fig6.Position=[10,5,17,8.5];
    ax.XLim=[0,max(x)];

    exportgraphics(fig6,[dirFigures,filesep,'statLevels',num2str(figidx),'.tiff']);


    %% Fluidization
    FG=timeseries2timetable(out.FG);        %Degree of fluidizations simulated
    gamma=timeseries2timetable(out.gamma);  %Mass diffusivity simulated


    %Set up figure
    fig7=figure(figidx+100);
    clf(fig7);
    ax=gca();
    box(ax,'on');


    %Excess fluidization
    plot(ax,x,FG{end,:}-1)
    ylabel(ax,'w_e/w_{mf} (-)');


    %Mass diffusivity
    hold(ax,'on');
    yyaxis(ax,'right');
    plot(ax,x,squeeze(gamma{end,:}))
    ylabel(ax,'Mass Diffusivity (m^2/s)');
    

    %Baffle positions
    xline(cumsum(xChambers(1:end-1)));
    hold(ax,'off');
    

    %Configure and save figure
    t7=title(ax,titleText,'Interpreter','latex');
    xlabel(ax,'Distance from Inlet (m)');

    fig7.Units='centimeters';
    fig7.Position=[10,5,17,8.5];
    ax.XLim=[0,max(x)];

    exportgraphics(fig7,[dirFigures,filesep,'statFluidization',num2str(figidx),'.tiff']);


    %% Export graphics for paper separately
    if figidx==21
        delete([t6,t7]);    %Delete titles

        %Resave as tiff
        exportgraphics(fig6,[dirFigures,filesep,'Figure6.tiff']);
        exportgraphics(fig7,[dirFigures,filesep,'Figure7.tiff']);

        %Save as eps
        exportgraphics(fig6,[dirFigures,filesep,'Figure6.eps']);
        exportgraphics(fig7,[dirFigures,filesep,'Figure7.eps']);
        
    end
end




