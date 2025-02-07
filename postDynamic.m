%% Post Processing of Dynamic Simulations
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
%Required products, version 24.1:
%   - MATLAB
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
    h4=dyn.h4+hBed;                     %Measured bed level h4


    %Residuals
    hcontInter=interp1(hsim.Time,hcont,dyn.Time);   %Simulated bed level at measured times
    resid=h4-hcontInter;
    
    
    %Set up figure
    fig6=figure(figidx);
    clf(fig6);
    ax=gca();
    box(ax,'on');
    hold(ax,'on');
    colors=ax.ColorOrder;

    plot(ax,dyn.Time,h4.*10^3);
    plot(ax,hsim.Time,hcont.*10^3);
    plot(ax,dyn.Time,dyn.AC1set.*10^3,'Color',colors(2,:),'LineStyle','--');
    hold(ax,'off');

    
    %Add residuals plot
    ax2=axes(fig6,'Units','centimeters','Position',[8,4,3,3]);
    histogram(ax2,resid.*10^3,'Normalization','pdf');

    hold(ax2,'on');

    x=linspace(ax2.XLim(1),ax2.XLim(2),100);

    pd=fitdist(resid.*10^3,'Normal');   %Probability distribution
    mu=pd.mu;                           %Expected value
    sigma=pd.sigma;                     %Standard deviation

    norm=exp(-(x-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));  %Normal distribution
    plot(x,norm,'LineWidth',1.5);

    hold(ax2,'off');

    [~,p]=chi2gof(resid);   %Chi-square goodness of fit p-value
    title(ax2,sprintf('\\mu=%.2f mm, p_\\chi=%.3f',mu,p));

    xlabel(ax2,'Residuals (mm)');
    ylabel(ax2,'PDF (-)');


    %Configure and save figure
    legend(ax,{'Measured','Predicted','Setpoint'});

    t6=title(ax,titleText,'Interpreter','latex');
    xlabel(ax,'Time (HH:MM:SS)');
    ylabel(ax,'Bed Level h_4 (mm)');

    ax.XLim=[0,max(dyn.Time)];
    ax.YLim=[ax.YLim(1),520];

    fig6.Units='centimeters';
    fig6.Position=[0.02,12.18,17,8.5];
    

    exportgraphics(fig6,[dirFigures,filesep,'stepRespContr',num2str(figidx),'.tiff'],...
        'Resolution',600);


    %% Other bed level responses
    hmeas=[dyn.h6,dyn.h5,dyn.h4];           %Measured bed levels
    hsimChambers=hsim{:,1}(:,posBedLevel);  %Predicted bed levels at each measurement position


    %Set up figure
    fig7=figure(figidx+100);
    clf(fig7);
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

    legend(ax,legItems,{'h_6','h_5','h_4','Measured','Predicted'},'Location','bestoutside');

    %Configure and save figure
    t7=title(ax,titleText,'Interpreter','latex');
    xlabel(ax,'Time (HH:MM:SS)');
    ylabel(ax,'Bed Level (m)');

    fig7.Units='centimeters';
    fig7.Position=[23.62,12.18,17,8.5];
    ax.XLim=[0,max(dyn.Time)];
    ax.YLim=[min(hmeas,[],'all')*0.85,max(hmeas,[],'all')*1.05]+hBed;

    ax2=axes(fig7,'Units','centimeters','Position',[12.5,1.25,3.6,3.59]);
    imshow([dirFigures,filesep,'StepResponseFigureInsert.tiff'],'Parent',ax2,'InitialMagnification','fit','Interpolation','bilinear','Reduce',false);

    exportgraphics(fig7,[dirFigures,filesep,'stepRespAll',num2str(figidx),'.tiff'],...
        'Resolution',600);


    %% Actuating value
    Ysim=timeseries2timetable(out.Ypid);
    Y=Ysim.Data(:,1);           %Predicted actuating value


    %Set up figure
    fig918=figure(figidx+900);
    clf(fig918);
    ax=gca();
    box(ax,'on');
    hold(ax,'on');

    plot(dyn.Time,dyn.AC1);
    plot(Ysim.Time,Y);

    %Mark start of step
    xline(ax,duration(0,0,20));
    hold(ax,'off');

    %Configure and save figure
    legend(ax,{'Measured','Predicted'});

    title(ax,titleText,'Interpreter','latex');
    xlabel(ax,'Time (HH:MM:SS)');
    ylabel(ax,'Valve Actuating Value (-)');

    fig918.Units='centimeters';
    fig918.Position=[0.02,0.83,17,8.5];
    ax.XLim=[0,max(dyn.Time)];

    exportgraphics(fig918,[dirFigures,filesep,'stepRespValve',num2str(figidx),'.tiff'],...
        'Resolution',600);


    %% Export graphics for paper separately
    if figidx==9
        delete([t6,t7]);    %Delete titles

        %Resave as tiff
        exportgraphics(fig6,[dirFigures,filesep,'Figure6.tiff'],'Resolution',600);
        exportgraphics(fig7,[dirFigures,filesep,'Figure7.tiff'],'Resolution',600);

        %Save as eps
        warning('off','MATLAB:print:ContentTypeImageSuggested');
        exportgraphics(fig6,[dirFigures,filesep,'Figure6.eps'],'ContentType','vector');
        exportgraphics(fig7,[dirFigures,filesep,'Figure7.eps']);
        warning('on','MATLAB:print:ContentTypeImageSuggested');

        %Save as fig
        savefig(fig6,[dirFigures,filesep,'Figure6']);
        savefig(fig7,[dirFigures,filesep,'Figure7']);
    end


end






