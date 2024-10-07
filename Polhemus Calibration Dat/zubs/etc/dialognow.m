

function dialognow(options)
FMSmax=66;
CPstr='CP'; onoff={'on','off'};
CPscell={'Control','Patient'};
LRscell={'Left','Right'};
LMacell={'Less','More'};

if nargin==0,
    options.Snum=[];
    options.CP=[];
    options.LR=[];
    options.LRdmg=[];
    options.LMa=[];
    options.FMS=[];
    options.position=[10 114 250 200]; %[-1833 114 250 150];
    options.done=0; end
if isempty(options.CP), options.CP=0; end
if isempty(options.FMS),
    if options.CP==1, options.FMS=nan; else options.FMS=50; end, end
if isempty(options.Snum), options.Snum=xtractSnum(options.CP); end
if isempty(options.LR), options.LR=1; end
if isempty(options.LRdmg), options.LRdmg=0; end
if isempty(options.LMa), options.LMa=0; end
if isempty(options.position), options.position=[300 300 250 150]; end

d=dialog('Position',options.position,'Name',['Sub: ' CPstr(options.CP+1) num2str(options.Snum)]);

%%%%%%%%%%%%%%%%%%%%
%%%%%% TEXT %%%%%%%%
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%%% C/P TEXT %%%
options.CPpos=[.3*options.position(3) options.position(4)-40 120 25];
CPtxt=uicontrol('Parent',d,...
    'Style','text',...
    'FontWeight','bold',...
    'HorizontalAlignment','left',...
    'FontSize',10,...
    'Position',options.CPpos,...
    'String',['SubType: ' CPscell{options.CP+1}]);

%%%%%%%%%%%%%%%%
%%% L/R TEXT %%%
options.LRpos=[.3*options.position(3) options.position(4)-60 120 25];
LRtxt=uicontrol('Parent',d,...
    'Style','text',...
    'FontWeight','bold',...
    'HorizontalAlignment','left',...
    'FontSize',10,...
    'Position',options.LRpos,...
    'String',['Handedness: ' LRscell{options.LR+1}]);

%%%%%%%%%%%%%%%%%%%%%%%
%%% HEM DAMAGE TEXT %%%
options.DMGpos=[.3*options.position(3) options.position(4)-80 120 25];
DMGtxt=uicontrol('Parent',d,...
    'Style','text',...
    'Visible',onoff{(~options.CP)+1},...
    'FontWeight','bold',...
    'HorizontalAlignment','left',...
    'FontSize',10,...
    'Position',options.DMGpos,...
    'String',['Hemisphere: ' LRscell{options.LRdmg+1}]);

%%%%%%%%%%%%%%%%%%%%%%%%
%%% LM AFFECTED TEXT %%%
options.LMApos=[.3*options.position(3) options.position(4)-100 120 25];
LMAtxt=uicontrol('Parent',d,...
    'Style','text',...
    'Visible',onoff{(~options.CP)+1},...
    'FontWeight','bold',...
    'HorizontalAlignment','left',...
    'FontSize',10,...
    'Position',options.LMApos,...
    'String',['Less/More: ' LMacell{options.LMa+1}]);

%%%%%%%%%%%%%%%%
%%% FMS TEXT %%%
options.FMSpos=[.3*options.position(3)-15 options.position(4)-125 55 25];
FMStxt=uicontrol('Parent',d,...
    'Style','text',...
    'Visible',onoff{(~options.CP)+1},...
    'FontWeight','bold',...
    'HorizontalAlignment','left',...
    'FontSize',10,...
    'Position',options.FMSpos,...
    'String',['FMS: ' num2str(options.FMS)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% UIcontrols %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
dH=15; dV=9;

%%%%%%%%%%%%%%%%%
%%% C/P CHECK %%%
uicontrol('Parent',d,...
    'pos',[options.CPpos(1)-dH options.CPpos(2)+dV 11 11],...
    'style','togglebutton',...
    'string','',...
    'value',options.CP, ...
    'background','white',...
    'foreground',[0 0 0], ...
    'Callback',@updateCP);

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HANDEDNESS TOGGLE %%%
uicontrol('Parent',d,...
    'pos',[options.LRpos(1)-dH options.LRpos(2)+dV 11 11],...
    'style','togglebutton',...
    'string','',...
    'value',options.LR, ...
    'background','white',...
    'foreground',[0 0 0], ...
    'Callback',@updateLR);

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HEM DAMAGE TOGGLE %%%
DMGb=uicontrol('Parent',d,...
    'style','togglebutton',...
    'Visible',onoff{(~options.CP)+1},...
    'pos',[options.DMGpos(1)-dH options.DMGpos(2)+dV 11 11],...
    'string','',...
    'value',options.LRdmg, ...
    'background','white',...
    'foreground',[0 0 0], ...
    'Callback',@updateDMG);

%%%%%%%%%%%%%%%%%%%%%%%%
%%% MORE/LESS TOGGLE %%%
LMAb=uicontrol('Parent',d,...
    'style','togglebutton',...
    'Visible',onoff{(~options.CP)+1},...
    'pos',[options.LMApos(1)-dH options.LMApos(2)+dV 11 11],...
    'string','',...
    'value',options.LMa, ...
    'background','white',...
    'foreground',[0 0 0], ...
    'Callback',@updateLMA);

%%%%%%%%%%%%%%%%%%
%%% FNS SLIDER %%%
FMSb=uicontrol('Parent',d,...
    'style','slider',...
    'Visible',onoff{(~options.CP)+1},...
    'pos',[options.position(3)/2-8 options.FMSpos(2)+4 80 20],...
    'Min',1,'Max',FMSmax,'Value',FMSmax-15,...
    'SliderStep',[1 5]/(FMSmax-1),...
    'Callback',@updateSLD);


%%%%%%%%%%%%%
%%% CLOSE %%%
uicontrol('Parent',d,...
    'String','done',...
    'Position',[options.position(3)/2-35 20 70 20],...
    'Callback',@alldone);




%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% CALLBACKS %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
    function alldone(~,~)
        save tmpsubdata options;
        retrievetempsubdata; 
        delete(gcf); end

    function updateSLD(source,~)
        options.FMS=source.Value;
        FMStxt.String=['FMS: ' num2str(options.FMS)]; end

    %%%% L/R HANDEDNESS
    function updateLR(source,~)
        options.LR=source.Value;
        LRtxt.String=['Handedness: ' LRscell{options.LR+1}]; end

    %%%% CONTROL / PATIENT   
    function updateCP(source,~)
        options.CP=source.Value;
        options.Snum=xtractSnum(options.CP);
        
        set(d,'Name',['Sub: ' CPstr(options.CP+1) num2str(options.Snum)]);
        CPtxt.String=['SubType: ' CPscell{options.CP+1}];
        DMGtxt.Visible=onoff{(~options.CP)+1}; DMGb.Visible=onoff{(~options.CP)+1}; 
        LMAtxt.Visible=onoff{(~options.CP)+1}; LMAb.Visible=onoff{(~options.CP)+1};
        FMStxt.Visible=onoff{(~options.CP)+1}; FMSb.Visible=onoff{(~options.CP)+1};  end


%%%% L/M AFFECTED
    function updateLMA(source,~)
        options.LMa=source.Value;
        LMAtxt.String=['L/M affected: ' LMacell{options.LMa+1}]; end

%%%% L/R HEMISPHERE DAMAGE
    function updateDMG(source,~)
        options.LRdmg=source.Value;
        DMGtxt.String=['Hemisphere: ' LRscell{options.LRdmg+1}]; end, end
