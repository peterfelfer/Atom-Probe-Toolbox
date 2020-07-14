function varargout = excessGUI(varargin)
% EXCESSGUI MATLAB code for excessGUI.fig
%      EXCESSGUI, by itself, creates a new EXCESSGUI or raises the existing
%      singleton*.
%
%      H = EXCESSGUI returns the handle to a new EXCESSGUI or the handle to
%      the existing singleton*.
%
%      EXCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXCESSGUI.M with the given input arguments.
%
%      EXCESSGUI('Property','Value',...) creates a new EXCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before excessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to excessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help excessGUI

% Last Modified by GUIDE v2.5 01-Nov-2013 14:54:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @excessGUI_OpeningFcn, ...
    'gui_OutputFcn',  @excessGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before excessGUI is made visible.
function excessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to excessGUI (see VARARGIN)

% Choose default command line output for excessGUI
handles.output = hObject;

%% variable setup
handles.cumulative = varargin{1};

sz = size(handles.cumulative);
if sz(2) < sz(1)
    handles.cumulative = handles.cumulative';
end

handles.interfaceArea = varargin{2};
handles.interfaceLoc = varargin{3};
numAtom = length(handles.cumulative);


% % downsampling if required
% MAXATOMS = 10000;
% 
% if numAtom > MAXATOMS
%     % downsample cumulative curve
%     handles.idx = round(linspace(1,numAtom,MAXATOMS));
%     handles.cumulative = handles.cumulative(handles.idx);
%     
% else
%     handles.idx = linspace(1,numAtom,numAtom);
% end
    



% if isSurf is true, only one approximation will be used
if (handles.interfaceLoc <=1) | (handles.interfaceLoc>=numAtom) % for surfaces/dislocations/clusters
    handles.isSurf = true;
else
    handles.isSurf = false;
end



%% plotting of cumulative curve
handles.cumulative_plothandle = plot(handles.cumulative,'Linewidth',2);
set(get(gca,'XLabel'),'String','cumulative number of atoms');
set(get(gca,'YLabel'),'String','cumulative number of species atoms');
hold all;



handles.x = 1:length(handles.cumulative);


%% intial limits

delta = round(numAtom/10);
loc = handles.interfaceLoc;

limL = loc-delta;
if limL<=1
    limL = loc;
elseif limL<loc/2
    limL = loc/2;
end

if handles.isSurf
    limL = 1;
end


limU = loc+delta;
if limU>numAtom
    limU=loc;
elseif (limU-loc)>((numAtom-loc)/2)
    limU = loc+(numAtom-loc)/2;
end


handles.limits = round([1 limL limU numAtom]);



%% plotting of limits and linear regression


handles.limits_plothandle = stem(handles.limits,repmat(max(handles.cumulative),[1 4]),...
    'k','Marker','none');


% linear regression (lower)
if ~handles.isSurf
    lin_reg_1 = polyfit(handles.x(handles.limits(1):handles.limits(2)),...
        handles.cumulative(handles.limits(1):handles.limits(2)),1);
    a1 = lin_reg_1(1);
    b1 = lin_reg_1(2);
else
    a1 = 0;
    b1 = 0;
end

hold on;
handles.loReg_plothandle = plot([handles.x(1)  handles.x(end)],[a1*handles.x(1)+b1 a1*handles.x(end)+b1],'k:','linewidth',1);


% linear regression (upper)
lin_reg_2 = polyfit(handles.x(handles.limits(3):handles.limits(4)),...
    handles.cumulative(handles.limits(3):handles.limits(4)),1);
a2 = lin_reg_2(1);
b2 = lin_reg_2(2);
hold on;
handles.upReg_plothandle = plot([handles.x(1)  handles.x(end)],[a2*handles.x(1)+b2 a2*handles.x(end)+b2],'k:','linewidth',1);



%% interfaceLoc
intersect1 = a1*handles.interfaceLoc+b1;
intersect2 = a2*handles.interfaceLoc+b2;

hold on;
handles.IEline_plothandle = plot([handles.interfaceLoc handles.interfaceLoc],[intersect1 intersect2],...
    '+--k');

%IE calculation
handles.IEcount = intersect2 - intersect1;
handles.IE = (intersect2 - intersect1) / handles.interfaceArea;


% partial excesses
intersectMax = handles.cumulative(round(handles.interfaceLoc));
% partial excess 1
handles.partialIE2 = (intersect2 - intersectMax) / handles.interfaceArea;
% partial excess 2
handles.partialIE1 = (intersectMax - intersect1) / handles.interfaceArea;


set(handles.excess_string,'String',['interfacial excess: ' num2str(handles.IE,3) ' at/nm2 (' num2str(round(handles.IEcount)) ' at)']);
set(handles.partialIE1_string,'String',['partial excess: ' num2str(handles.partialIE1,3) ' at/nm2']);
set(handles.partialIE2_string,'String',['partial excess: ' num2str(handles.partialIE2,3) ' at/nm2']);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes excessGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = excessGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in changeLimits_button.
function changeLimits_button_Callback(hObject, eventdata, handles)
% hObject    handle to changeLimits_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% get graphical input for new limit

[x y] = ginput(2);

dist = abs(handles.limits - x(1));
target = find(dist == min(dist));
target = target(1); %in case you click right in the middle between two limits

if handles.isSurf & ((target == 1) | (target == 2))
    target = 3;
end

handles.limits(target) = round(x(2));

recalculateValues(hObject,handles);





function recalculateValues(hObject,handles)
set(handles.limits_plothandle,'XData',handles.limits);

% linear regression (lower)
if ~handles.isSurf
    lin_reg_1 = polyfit(handles.x(handles.limits(1):handles.limits(2)),...
        handles.cumulative(handles.limits(1):handles.limits(2)),1);
    a1 = lin_reg_1(1);
    b1 = lin_reg_1(2);
else
    a1 = 0;
    b1 = 0;
end


hold on;
set(handles.loReg_plothandle,'YData',[a1*handles.x(1)+b1 a1*handles.x(end)+b1])

% linear regression (upper)
lin_reg_2 = polyfit(handles.x(handles.limits(3):handles.limits(4)),...
    handles.cumulative(handles.limits(3):handles.limits(4)),1);
a2 = lin_reg_2(1);
b2 = lin_reg_2(2);
hold on;
set(handles.upReg_plothandle,'YData',[a2*handles.x(1)+b2 a2*handles.x(end)+b2]);


%% interfaceLoc
intersect1 = a1*handles.interfaceLoc+b1;
intersect2 = a2*handles.interfaceLoc+b2;

hold on;
set(handles.IEline_plothandle,'YData',[intersect1 intersect2]);
set(handles.IEline_plothandle,'XData',[handles.interfaceLoc handles.interfaceLoc]);

%IE calculation
handles.IEcount = intersect2 - intersect1;
handles.IE = (intersect2 - intersect1) / handles.interfaceArea;


% partial excesses
intersectMax = handles.cumulative(round(handles.interfaceLoc));
% partial excess 1
handles.partialIE2 = (intersect2 - intersectMax) / handles.interfaceArea;
% partial excess 2
handles.partialIE1 = (intersectMax - intersect1) / handles.interfaceArea;


set(handles.excess_string,'String',['interfacial excess: ' num2str(handles.IE,3) ' at/nm2 (' num2str(round(handles.IEcount)) ' at)']);
set(handles.partialIE1_string,'String',['partial excess: ' num2str(handles.partialIE1,3) ' at/nm2']);
set(handles.partialIE2_string,'String',['partial excess: ' num2str(handles.partialIE2,3) ' at/nm2']);

guidata(hObject,handles)


% --- Executes on button press in interfaceMax_button.
function interfaceMax_button_Callback(hObject, eventdata, handles)
% hObject    handle to interfaceMax_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of interfaceMax_button


% --------------------------------------------------------------------
function changeLimits_Callback(hObject, eventdata, handles)
% hObject    handle to changeLimits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
changeLimits_button_Callback(hObject, eventdata, handles)

guidata(hObject,handles);


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


switch eventdata.Key
    
    case 'l'
        %changeLimits_button_Callback(hObject,eventdata,handles)
        
end



guidata(hObject,handles);


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uiputfile('*.fig','save figure');

plotHandles = [handles.cumulative_plothandle handles.limits_plothandle...
    handles.loReg_plothandle handles.upReg_plothandle...
    handles.IEline_plothandle];


ftmp = figure;
set(gcf,'Name','cumulative diagram');
set(gcf,'Color',[1 1 1]);

atmp = axes;
copyobj(plotHandles, atmp);


txt1 = ['interfacial excess: ' num2str(handles.IE,3) ' at/nm2 (' num2str(round(handles.IEcount)) ' at)'];

text(length(handles.cumulative)*0.1,max(handles.cumulative)*0.95,txt1);
txt2 = ['partials: IE1 =' num2str(handles.partialIE1,3) ' at/nm2'...
    ',  IE2 = ' num2str(handles.partialIE2,3) ' at/nm2'];
text(length(handles.cumulative)*0.1,max(handles.cumulative)*0.8,txt2);


set(get(gca,'XLabel'),'String','cumulative number of atoms');
set(get(gca,'YLabel'),'String','cumulative number of excess species');

if file
    saveas(ftmp, [path file]);
    delete(ftmp);
    
else
    plotedit on
    
end
