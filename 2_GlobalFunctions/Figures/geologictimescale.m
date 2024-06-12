function geologictimescale(xmin,xmax,Yaxisopt,Xaxisopt,axisselect,...
                            colorselect,divisions,labels,sf,axnum,fontname,...
                            fontsize,reccolor,textcolor)
                        
                        
%SCRIPT TO ADD GEOLOGIC TIME SCALE TO FIGURES
%Created: 08/2020 (E. Judd)
%Last updated: 08/2023 (E. Judd)

%INPUTS:
%xmin: minimum age of plot (in myrs)
    %default: 0 Ma
%xmax: maximum age of plot (in myrs)
    %default: 538.8 Ma
%Yaxisopts: do you want the y-axis normal or reverse (L-R)?
    %'normal' (default)
    %'reverse'
%Xaxisopts: do you want the x-axis normal or reverse?
    %'normal' (defualt; time proceeds right to left)
    %'reverse' (time proceeds left to right)
%axisselect: which axis do you want to plot the time scale on?
%            (default: gca)
%colorselect: do you want to use the standard timescale colors or "fun"
%             colors? 
    %'standard' (default)
    %'fun' (currently in beta mode; fun only works in Cenozoic)
%divisions: the finest resolution age division you'd like to show
    %'periods' (default, plots Eras and Periods)
    %'epochs' (plots Periods and Epochs)
    %'stages' (plots Periods and Stages)
%labels: should the finest resultion age division include labels?
    %'on' (default)
    %'off' (useful if you're plotting stages)
%sf: scale factor: scalar used to adjust the width of the time scale
    % scalar (default: 7.5)
%axnum: number of axes 
    %1 (default; ticks added to right y-axis)
    %2 (no ticks added to right y-axis)
%fontname
    %'helvetica' (default)
    %any other accepted matlab font
%fontsize
    %scalar (default: 11)

%SYNTAX:
%geologictimescale            <--- uses all defaults
%geologictimescale(xmin,xmax) <--- fixes xlim otherwise uses defaults
%geologictimescale(...)       <--- user specifies all inputs
    
    
% (1) Load GTS file (CHANGE BASED ON FILENAME/LOCATION)
    load("GTS2020.mat")
% (2) Specify defaults basaed on nargin
    if nargin == 0
        xmin = GTS.UpperBoundary(1);
        xmax = GTS.LowerBoundary(end);
        Yaxisopt = 'normal';
        Xaxisopt = 'normal';
        axisselect = gca;
        colorselect = 'standard';
        divisions = 'periods';
        labels = 'on';
        sf = 7.5;
        axnum = 1;
        fontname = 'helvetica';
        fontsize = 11;
        reccolor = [0 0 0];
        textcolor = [0 0 0];
    elseif nargin == 2
        Yaxisopt = 'normal';
        Xaxisopt = 'normal';
        axisselect = gca;
        colorselect = 'standard';
        divisions = 'periods';
        labels = 'on';
        sf = 7.5;
        axnum = 1;
        fontname = 'helvetica';
        fontsize = 11;
        reccolor = [0 0 0];
        textcolor = [0 0 0];
    elseif nargin == 3
        Xaxisopt = 'normal';
        axisselect = gca;
        colorselect = 'standard';
        divisions = 'periods';
        labels = 'on';
        sf = 7.5;
        axnum = 1;
        fontname = 'helvetica';
        fontsize = 11;
        reccolor = [0 0 0];
        textcolor = [0 0 0];
    elseif nargin == 4
        axisselect = gca;
        colorselect = 'standard';
        divisions = 'periods';
        labels = 'on';
        sf = 7.5;
        axnum = 1;
        fontname = 'helvetica';
        fontsize = 11;
        reccolor = [0 0 0];
        textcolor = [0 0 0];
    elseif nargin == 5
        colorselect = 'standard';
        divisions = 'periods';
        labels = 'on';
        sf = 7.5;
        axnum = 1;
        fontname = 'helvetica';
        fontsize = 11;
        reccolor = [0 0 0];
        textcolor = [0 0 0];
    elseif nargin == 6
        divisions = 'periods';
        labels = 'on';
        sf = 7.5;
        axnum = 1;
        fontname = 'helvetica';
        fontsize = 11;
        reccolor = [0 0 0];
        textcolor = [0 0 0];
    elseif nargin == 7
        labels = 'on';
        sf = 7.5;
        axnum = 1;
        fontname = 'helvetica';
        fontsize = 11;
        reccolor = [0 0 0];
        textcolor = [0 0 0];
    elseif nargin == 8
        sf = 7.5;
        axnum = 1;
        fontname = 'helvetica';
        fontsize = 11;
        reccolor = [0 0 0];
        textcolor = [0 0 0];
    elseif nargin == 9
        axnum = 1;
        fontname = 'helvetica';
        fontsize = 11;
        reccolor = [0 0 0];
        textcolor = [0 0 0];
    elseif nargin == 10
        fontname = 'helvetica';
        fontsize = 11;
        reccolor = [0 0 0];
        textcolor = [0 0 0];
    elseif nargin == 10
        fontsize = 11;
        reccolor = [0 0 0];
        textcolor = [0 0 0];
    elseif nargin == 11
        reccolor = [0 0 0];
        textcolor = [0 0 0];
    elseif nargin == 12
        textcolor = [0 0 0];
    end
% (3) Define colors, times, and names of age divisions
    % (a) Eras:
    Color.Eras = [242,249,29;103,197,202;153,192,141]/255;
    Time.Eras = [0.0;66.04;251.9;538.8];
    Name.Eras = ["Cenozoic";"Mesozoic";"Paleozoic"];
    % (b) Periods:
    Color.Periods = [GTS.Rperiod,GTS.Gperiod,GTS.Bperiod];
    [~, idx] = unique(Color.Periods,'rows');idx = sort(idx);
    Color.Periods = Color.Periods(idx,:);
    if strcmpi(colorselect,'fun') == 1
    Color.Periods = [hex2rgb('#868e71',1);hex2rgb('#899099',1);hex2rgb('#8c7a84',1);127,198,78;52,178,201;129,43,146;...
        240,64,40;123,165,153;203,140,55;179,225,182;0,146,112;127,160,86];
    end    
    Time.Periods = [GTS.UpperBoundary(idx);GTS.LowerBoundary(end)];
    Name.Periods = GTS.Period(idx);Name.Periods(1) = {'Q'};
    Abbrev.Periods = ["Q";"N";"Pg";"K";"J";"T";"P";"C";"D";"S";"O";"Cm"];
    % (c) Epochs:
    Color.Epochs = [GTS.Repoch,GTS.Gepoch,GTS.Bepoch];
    [~, idx] = unique(Color.Epochs,'rows');idx = sort(idx);
    Color.Epochs = Color.Epochs(idx,:);
    if strcmpi(colorselect,'fun')
        Color.Epochs = hex2rgb({'#9ba58c','#717856', '#acb5bb','#656b77', '#dac3b8','#a58e91','#736678'},1);
    end
    Time.Epochs = [GTS.UpperBoundary(idx);GTS.LowerBoundary(end)];
    Name.Epochs = GTS.Epoch(idx);Name.Epochs(1) = {'H'};
        Name.Epochs(2) = {'Pleist'};Name.Epochs(3) = {'Plio'};
    % (d) Stages:
    Color.Stages = [GTS.Rstage,GTS.Gstage,GTS.Bstage];
    Time.Stages = [GTS.UpperBoundary;GTS.LowerBoundary(end)];
    Name.Stages = GTS.Stage;

%(4) Adjust range of starting parameters to data
    %(a) Eras
    idx1 = find(Time.Eras<=xmin);idx2=find(Time.Eras>=xmax);
    Time.Eras=Time.Eras(idx1(end):idx2(1));Time.Eras(1)=xmin;Time.Eras(end)=xmax;
    idx2=idx2(1)-1;
    Name.Eras=Name.Eras(idx1(end):idx2(1));
    Color.Eras=Color.Eras(idx1(end):idx2(1),:);
    %(b) Periods
    idx1 = find(Time.Periods<=xmin);idx2=find(Time.Periods>=xmax);
    Time.Periods=Time.Periods(idx1(end):idx2(1));Time.Periods(1)=xmin;Time.Periods(end)=xmax;
    idx2=idx2(1)-1;
    Name.Periods=Name.Periods(idx1(end):idx2(1));
    Color.Periods=Color.Periods(idx1(end):idx2(1),:);
    Abbrev.Periods=Abbrev.Periods(idx1(end):idx2(1));
    %(c) Epochs
    idx1 = find(Time.Epochs<=xmin);idx2=find(Time.Epochs>=xmax);
    Time.Epochs=Time.Epochs(idx1(end):idx2(1));Time.Epochs(1)=xmin;Time.Epochs(end)=xmax;
    idx2=idx2(1)-1;
    Name.Epochs=Name.Epochs(idx1(end):idx2(1));
    Color.Epochs=Color.Epochs(idx1(end):idx2(1),:);
    %(d) Stages
    idx1 = find(Time.Stages<=xmin);idx2=find(Time.Stages>=xmax);
    Time.Stages=Time.Stages(idx1(end):idx2(1));Time.Stages(1)=xmin;Time.Stages(end)=xmax;
    idx2=idx2(1)-1;
    Name.Stages=Name.Stages(idx1(end):idx2(1));
    Color.Stages=Color.Stages(idx1(end):idx2(1),:);
%(5) Define starting paramets for figure
    hold on,ax = axisselect;yl = ax.YLim;yt = ax.YTick;xt = ax.XTick;
    ax.XLim = [xmin,xmax];
    if strcmpi(Yaxisopt,'normal')
        ax.YLim = [yl(1)-range(yl)/sf yl(2)];
    elseif strcmpi(Yaxisopt,'reverse')
        set(ax, 'Ydir', 'reverse')
        yyaxis right, set(ax, 'Ydir', 'reverse'), yyaxis left
        ax.YLim = [yl(1) yl(2)+range(yl)/sf];
    end
    ax.YTick = yt;
    if strcmpi(Xaxisopt,'reverse')
        set(ax, 'Xdir', 'reverse')
    end
    ym = ax.YLim;
    
    if divisions == "all"
        recheight = range(yl)/(sf*3);
    else
        recheight = range(yl)/(sf*2);
    end
    
%(6) Add rectangles and text for each division
    %(Option 1): If plotting by periods
    %(a) Add Era lables
    if strcmpi(divisions,'periods')==1
    for ii = 1:length(Name.Eras)
    if strcmpi(Yaxisopt,'normal')
        rectangle(ax,'Position',[Time.Eras(ii),ym(1),range(Time.Eras(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Eras(ii,:));
        text(ax,Time.Eras(ii)+range(Time.Eras(ii:ii+1))/2,ym(1)+recheight/2,Name.Eras(ii),...
            'HorizontalAlignment','center','FontName',fontname,'color',textcolor)
    elseif strcmpi(Yaxisopt,'reverse')
        rectangle(ax,'Position',[Time.Eras(ii),yl(2)+range(yl)/(sf*2),range(Time.Eras(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Eras(ii,:));
        text(ax,Time.Eras(ii)+range(Time.Eras(ii:ii+1))/2,ym(1)+recheight/2,Name.Eras(ii),...
            'HorizontalAlignment','center','FontName',fontname,'color',textcolor,'VerticalAlignment','middle')
    end
    end    
    %(b) Add Period labels
    for ii = 1:length(Name.Periods)
    if strcmpi(Yaxisopt,'normal')
        rectangle(ax,'Position',[Time.Periods(ii),yl(1)-range(yl)/sf+range(yl)/(sf*2),...
            range(Time.Periods(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Periods(ii,:));
        if strcmpi(labels,'on')
            text(ax,Time.Periods(ii)+range(Time.Periods(ii:ii+1))/2,yl(1)-range(yl)/sf+range(yl)/(sf*2)+recheight/2,...
            Name.Periods(ii),...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontName',fontname,'color',textcolor)
        end
    elseif strcmpi(Yaxisopt,'reverse')
        rectangle(ax,'Position',[Time.Periods(ii),yl(2),...
            range(Time.Periods(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Periods(ii,:));
        if strcmpi(labels,'on')
            text(ax,Time.Periods(ii)+range(Time.Periods(ii:ii+1))/2,recheight+recheight/2,Name.Periods(ii),...
            'HorizontalAlignment','center','FontName',fontname,'color',textcolor,'VerticalAlignment','middle')
        end
    end
    end 

    elseif strcmpi(divisions,'epochs')==1
    %(Option 2): If plotting by epoch
    %(a) Add Period lables
    for ii = 1:length(Name.Periods)
    if strcmpi(Yaxisopt,'normal')
        rectangle(ax,'Position',[Time.Periods(ii),ym(1),range(Time.Periods(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Periods(ii,:));
        text(ax,Time.Periods(ii)+range(Time.Periods(ii:ii+1))/2,ym(1)+recheight/2,Abbrev.Periods(ii),...
            'HorizontalAlignment','center','FontName',fontname,'color',textcolor)
    elseif strcmpi(Yaxisopt,'reverse')
        rectangle(ax,'Position',[Time.Periods(ii),yl(2)+range(yl)/(sf*2),range(Time.Periods(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Periods(ii,:));
        text(ax,Time.Periods(ii)+range(Time.Periods(ii:ii+1))/2,ym(1)+recheight/2,Abbrev.Periods(ii),...
            'HorizontalAlignment','center','FontName',fontname,'color',textcolor,'VerticalAlignment','middle')
    end
    end    
    %(b) Add Epoch labels
    for ii = 1:length(Name.Epochs)
    if strcmpi(Yaxisopt,'normal')
        rectangle(ax,'Position',[Time.Epochs(ii),yl(1)-range(yl)/sf+range(yl)/(sf*2),...
            range(Time.Epochs(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Epochs(ii,:));
        if strcmpi(labels,'on')
            text(ax,Time.Epochs(ii)+range(Time.Epochs(ii:ii+1))/2,yl(1)-range(yl)/sf+range(yl)/(sf*2)+recheight/2,...
            Name.Epochs(ii),...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontName',fontname,'color',textcolor)
        end
    elseif strcmpi(Yaxisopt,'reverse')
        rectangle(ax,'Position',[Time.Epochs(ii),yl(2),...
            range(Time.Epochs(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Epochs(ii,:));
        if strcmpi(labels,'on')
            text(ax,Time.Epochs(ii)+range(Time.Epochs(ii:ii+1))/2,recheight+recheight/2,Name.Epochs(ii),...
            'HorizontalAlignment','center','FontName',fontname,'color',textcolor,'VerticalAlignment','middle')
        end
    end
    end
    
    elseif strcmpi(divisions,'stages')==1
    %(Option 3): If plotting by stage
    %(a) Add Period lables
    for ii = 1:length(Name.Periods)
    if strcmpi(Yaxisopt,'normal')
        rectangle(ax,'Position',[Time.Periods(ii),ym(1),range(Time.Periods(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Periods(ii,:));
        text(ax,Time.Periods(ii)+range(Time.Periods(ii:ii+1))/2,ym(1)+recheight/2,Abbrev.Periods(ii),...
            'HorizontalAlignment','center','FontName',fontname,'color',textcolor)
    elseif strcmpi(Yaxisopt,'reverse')
        rectangle(ax,'Position',[Time.Periods(ii),yl(2)+range(yl)/(sf*2),range(Time.Periods(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Periods(ii,:));
        text(ax,Time.Periods(ii)+range(Time.Periods(ii:ii+1))/2,ym(1)+recheight/2,Abbrev.Periods(ii),...
            'HorizontalAlignment','center','FontName',fontname,'color',textcolor,'VerticalAlignment','middle')
    end
    end    
    %(b) Add Stage labels
    for ii = 1:length(Name.Stages)
    if strcmpi(Yaxisopt,'normal')
        rectangle(ax,'Position',[Time.Stages(ii),yl(1)-range(yl)/sf+range(yl)/(sf*2),...
            range(Time.Stages(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Stages(ii,:));
        if strcmpi(labels,'on')
            lab = char(Name.Stages(ii));
            text(ax,Time.Stages(ii)+range(Time.Stages(ii:ii+1))/2,yl(1)-range(yl)/sf+range(yl)/(sf*2)+recheight/2,...
            lab(1:3),...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontName',fontname,'color',textcolor)
        end
    elseif strcmpi(Yaxisopt,'reverse')
        rectangle(ax,'Position',[Time.Stages(ii),yl(2),...
            range(Time.Stages(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Stages(ii,:));
        if strcmpi(labels,'on')
            text(ax,Time.Stages(ii)+range(Time.Stages(ii:ii+1))/2,recheight+recheight/2,Name.Stages(ii),...
            'HorizontalAlignment','center','FontName',fontname,'color',textcolor,'VerticalAlignment','middle')
        end
    end
    end 
    
    elseif strcmpi(divisions,'all')==1
    %(Option 4): If plotting by stage, period, and era (i.e., "all")
    %(a) Add era lables
    for ii = 1:length(Name.Eras)
    if strcmpi(Yaxisopt,'normal')
        rectangle(ax,'Position',[Time.Eras(ii),ym(1),range(Time.Eras(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Eras(ii,:));
        text(ax,Time.Eras(ii)+range(Time.Eras(ii:ii+1))/2,ym(1)+recheight/2,Name.Eras(ii),...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontName',fontname,'color',textcolor)
    elseif strcmpi(Yaxisopt,'reverse')
        rectangle(ax,'Position',[Time.Eras(ii),yl(2)+range(yl)/(sf*2),range(Time.Eras(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Eras(ii,:));
        text(ax,Time.Eras(ii)+range(Time.Eras(ii:ii+1))/2,ym(1)+recheight/2,Name.Eras(ii),...
            'HorizontalAlignment','center','FontName',fontname,'color',textcolor,'VerticalAlignment','middle')
    end
    end    
    %(b) Add Period lables
    for ii = 1:length(Name.Periods)
    if strcmpi(Yaxisopt,'normal')
        rectangle(ax,'Position',[Time.Periods(ii),ym(1)+recheight,range(Time.Periods(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Periods(ii,:));
        text(ax,Time.Periods(ii)+range(Time.Periods(ii:ii+1))/2,ym(1)+recheight+recheight/2,Abbrev.Periods(ii),...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontName',fontname,'color',textcolor)
    elseif strcmpi(Yaxisopt,'reverse')
        rectangle(ax,'Position',[Time.Periods(ii),yl(2)+range(yl)/(sf*2),range(Time.Periods(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Periods(ii,:));
        text(ax,Time.Periods(ii)+range(Time.Periods(ii:ii+1))/2,ym(1)+recheight/2,Abbrev.Periods(ii),...
            'HorizontalAlignment','center','FontName',fontname,'color',textcolor,'VerticalAlignment','middle')
    end
    end    
    %(b) Add Stage labels
    for ii = 1:length(Name.Stages)
    if strcmpi(Yaxisopt,'normal')
        rectangle(ax,'Position',[Time.Stages(ii),ym(1)+2*recheight,...
            range(Time.Stages(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Stages(ii,:));
        if strcmpi(labels,'on')
            lab = char(Name.Stages(ii));
            text(ax,Time.Stages(ii)+range(Time.Stages(ii:ii+1))/2,yl(1)-range(yl)/sf+range(yl)/(sf*2)+recheight/2,...
            lab(1:3),...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontName',fontname,'color',textcolor)
        end
    elseif strcmpi(Yaxisopt,'reverse')
        rectangle(ax,'Position',[Time.Stages(ii),yl(2),...
            range(Time.Stages(ii:ii+1)),recheight],...
            'EdgeColor',reccolor,'FaceColor',Color.Stages(ii,:));
        if strcmpi(labels,'on')
            text(ax,Time.Stages(ii)+range(Time.Stages(ii:ii+1))/2,recheight+recheight/2,Name.Stages(ii),...
            'HorizontalAlignment','center','FontName',fontname,'color',textcolor,'VerticalAlignment','middle')
        end
    end
    end   
    end
    
    set(findall(gca,'-property','FontSize'),'FontSize',fontsize)

%(7) Make cosmetic adjustments (box around plot, custom tick marks)
    box off
    colororder({'k','k'})
    yl=ax.YLim;xl=ax.XLim;
    rectangle(ax,'Position',[xl(1),yl(1),range(xl),range(yl)],'EdgeColor',reccolor,'LineWidth',1);
    ax.XAxis.TickDirection='out';ax.XAxis.MinorTick='on';
    if axnum == 1
        yt = ax.YTick;
        ytm = ax.YAxis.MinorTickValues;
        yyaxis(ax,'right');
        ax.YLim=yl;
        ax.YTick=yt;
        ax.YAxis(2).MinorTickValues=ytm;
        ax.YAxis(2).MinorTick=ax.YAxis(1).MinorTick;
        ax.YColor = 'k';
        ax.YTickLabel=[];yyaxis left
    end
end
