%% Post Processing of Dynamic Simulations
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
%All required files for this function can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.7948224
%
%
%
%This function anaylzes the results of the dynamic simulations conducted
%with the script "calcDynamic" and creates all published figures.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%Necessary files, classes, functions, and scripts:
%   - StepResponseFigureInsert.tiff in the dirFigures folder


function postDynamic(out,flow,dyn,hBed,posBedLevel,figidx,dirFigures)
    %% Common parameters
    %Figure title:
    titleText=['Step Response, Test ',num2str(figidx),...
                ', $w_{e}$/$w_{mf}$=',num2str(round(flow.FG2-1,1))];

    
    %% Controller response
    hsim=timeseries2timetable(out.h);
    hcont=hsim{:,1}(:,posBedLevel(3));  %Controlled (simulated) bed level h_4
    
    
    %Set up figure
    fig9=figure(figidx);
    clf(fig9);
    ax=gca();
    box(ax,'on');
    hold(ax,'on');
    colors=ax.ColorOrder;

    plot(ax,dyn.Time,dyn.h4+hBed);
    plot(ax,hsim.Time,hcont);
    plot(ax,dyn.Time,dyn.AC1set,'Color',colors(2,:),'LineStyle','--');
    hold(ax,'off');

    %Configure and save figure
    legend(ax,{'Measured','Simulated','Setpoint'});

    t9=title(ax,titleText,'Interpreter','latex');
    xlabel(ax,'Time (HH:MM:SS)');
    ylabel(ax,'Bed Level h_4 (m)');

    fig9.Units='centimeters';
    fig9.Position=[0.02,12.18,17,8.5];
    ax.XLim=[0,max(dyn.Time)];

    exportgraphics(fig9,[dirFigures,filesep,'stepRespContr',num2str(figidx),'.tiff']);


    %% Actuating value
    Ysim=timeseries2timetable(out.Ypid);
    Y=Ysim.Data(:,1);           %Simulated actuating value


    %Set up figure
    fig10=figure(figidx+100);
    clf(fig10);
    ax=gca();
    box(ax,'on');
    hold(ax,'on');

    plot(dyn.Time,dyn.AC1);
    plot(Ysim.Time,Y);

    %Mark start of step
    xline(ax,duration(0,0,20));
    hold(ax,'off');

    %Configure and save figure
    legend(ax,{'Measured','Simulated'});

    t10=title(ax,titleText,'Interpreter','latex');
    xlabel(ax,'Time (HH:MM:SS)');
    ylabel(ax,'Valve Actuating Value (-)');

    fig10.Units='centimeters';
    fig10.Position=[0.02,0.83,17,8.5];
    ax.XLim=[0,max(dyn.Time)];

    exportgraphics(fig10,[dirFigures,filesep,'stepRespValve',num2str(figidx),'.tiff']);


    %% Other bed level responses
    hmeas=[dyn.h6,dyn.h5,dyn.h4];           %Measured bed levels
    hsimChambers=hsim{:,1}(:,posBedLevel);  %Simulated bed levels at each measurement position


    %Set up figure
    fig11=figure(figidx+200);
    clf(fig11);
    ax=gca();
    box(ax,'on');
    hold(ax,'on');
    colors=ax.ColorOrder;

    for i=1:size(hmeas,2)
        plot(ax,dyn.Time,hmeas(:,i)+hBed,'Color',colors(i,:),'LineStyle','-');
        plot(ax,hsim.Time,hsimChambers(:,i),'Color',colors(i,:),'LineStyle','--');
    end

    %Mark start of step
    xline(ax,duration(0,0,20));

    %Create legend items
    legItems=repmat(line(),5,1);
    d0=duration(NaN,NaN,NaN);
    legItems(1)=plot(ax,d0,0,'Color',colors(1,:),'LineStyle','-');
    legItems(2)=plot(ax,d0,0,'Color',colors(2,:),'LineStyle','-');
    legItems(3)=plot(ax,d0,0,'Color',colors(3,:),'LineStyle','-');
    legItems(4)=plot(ax,d0,0,'Color','k','LineStyle','-');
    legItems(5)=plot(ax,d0,0,'Color','k','LineStyle','--');
    hold(ax,'off');

    legend(ax,legItems,{'h_6','h_5','h_4','Measured','Simulated'},'Location','bestoutside');

    %Configure and save figure
    t11=title(ax,titleText,'Interpreter','latex');
    xlabel(ax,'Time (HH:MM:SS)');
    ylabel(ax,'Bed Level (m)');

    fig11.Units='centimeters';
    fig11.Position=[23.62,12.18,17,8.5];
    ax.XLim=[0,max(dyn.Time)];
    ax.YLim=[min(hmeas,[],'all')*0.85,max(hmeas,[],'all')*1.05]+hBed;

    ax2=axes(fig11,'Units','centimeters','Position',[12.5,1.25,3.6,3.59]);
    imshow([dirFigures,filesep,'StepResponseFigureInsert.tiff'],'Parent',ax2,'InitialMagnification','fit','Interpolation','bilinear','Reduce',false);

    exportgraphics(fig11,[dirFigures,filesep,'stepRespAll',num2str(figidx),'.tiff']);


    %% Export graphics for paper separately
    if figidx==9
        delete([t9,t10,t11]);    %Delete titles

        %Resave as tiff
        exportgraphics(fig9,[dirFigures,filesep,'Figure9.tiff']);
        exportgraphics(fig10,[dirFigures,filesep,'Figure10.tiff']);
        exportgraphics(fig11,[dirFigures,filesep,'Figure11.tiff']);

        %Save as eps
        exportgraphics(fig9,[dirFigures,filesep,'Figure9.eps']);
        exportgraphics(fig10,[dirFigures,filesep,'Figure10.eps']);
        exportgraphics(fig11,[dirFigures,filesep,'Figure11.eps']);
    end


end






