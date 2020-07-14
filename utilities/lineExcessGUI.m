function varargout = lineExcessGUI(varargin)
% LINEEXCESSGUI MATLAB code for lineExcessGUI.fig
%      LINEEXCESSGUI, by itself, creates a new LINEEXCESSGUI or raises the existing
%      singleton*.
%
%      H = LINEEXCESSGUI returns the handle to a new LINEEXCESSGUI or the handle to
%      the existing singleton*.
%
%      LINEEXCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LINEEXCESSGUI.M with the given input arguments.
%
%      LINEEXCESSGUI('Property','Value',...) creates a new LINEEXCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lineExcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lineExcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lineExcessGUI

% Last Modified by GUIDE v2.5 21-Nov-2017 09:00:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lineExcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @lineExcessGUI_OutputFcn, ...
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


% --- Executes just before lineExcessGUI is made visible.
function lineExcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lineExcessGUI (see VARARGIN)

% Choose default command line output for lineExcessGUI
handles.output = hObject;


%% variable setup
handles.cumulative = varargin{1};

sz = size(handles.cumulative);
if sz(2) < sz(1)
    handles.cumulative = handles.cumulative';
end

handles.lineLength = varargin{2};
numAtom = length(handles.cumulative);


%% plotting of cumulative curve
handles.cumulative_plothandle = plot(handles.cumulative,'Linewidth',2);
set(get(gca,'XLabel'),'String','cumulative number of atoms');
set(get(gca,'YLabel'),'String','cumulative number of species atoms');
hold all;

handles.x = 1:length(handles.cumulative);



%% intial limits
limU = round(numAtom/10);

if limU > numAtom
    limU=loc;
elseif limU > (numAtom/2)
    limU = numAtom/2;
end


handles.limits = round([limU numAtom]);


%% plotting of limits and linear regression
handles.limits_plothandle = stem(handles.limits,repmat(max(handles.cumulative),[1 2]),...
    'k','Marker','none');


% linear regression
lin_reg = polyfit(handles.x(handles.limits(1):handles.limits(2)),...
    handles.cumulative(handles.limits(1):handles.limits(2)),1);
a = lin_reg(1);
b = lin_reg(2);
hold on;
handles.upReg_plothandle = plot([handles.x(1)  handles.x(end)],[a*handles.x(1)+b a*handles.x(end)+b],'k:','linewidth',1);

%% interfaceLoc
handles.IEcount = b;
handles.IE = b / handles.lineLength;

set(handles.excess_string,'String',['line excess: ' num2str(handles.IE,3) ' at/nm (' num2str(round(handles.IEcount)) ' at)']);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lineExcessGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);





% --- Outputs from this function are returned to the command line.
function varargout = lineExcessGUI_OutputFcn(hObject, eventdata, handles) 
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

handles.limits(target) = round(x(2));


recalculateValues(hObject,handles);



function recalculateValues(hObject,handles)
set(handles.limits_plothandle,'XData',handles.limits);


% linear regression (upper)
lin_reg = polyfit(handles.x(handles.limits(1):handles.limits(2)),...
    handles.cumulative(handles.limits(1):handles.limits(2)),1);
a = lin_reg(1);
b = lin_reg(2);
hold on;
set(handles.upReg_plothandle,'YData',[a*handles.x(1)+b a*handles.x(end)+b]);


%% interfaceLoc
intersect = b;

hold on;

%IE calculation
handles.IEcount = intersect;
handles.IE = intersect / handles.lineLength;

set(handles.excess_string,'String',['line excess: ' num2str(handles.IE,3) ' at/nm (' num2str(round(handles.IEcount)) ' at)']);

guidata(hObject,handles)




% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uiputfile('*.fig','save figure');

plotHandles = [handles.cumulative_plothandle handles.limits_plothandle...
    handles.upReg_plothandle];


ftmp = figure;
set(gcf,'Name','cumulative diagram');
set(gcf,'Color',[1 1 1]);

atmp = axes;
copyobj(plotHandles, atmp);


txt1 = ['line excess: ' num2str(handles.IE,3) ' at/nm (' num2str(round(handles.IEcount)) ' at)'];

text(length(handles.cumulative)*0.1,max(handles.cumulative)*0.95,txt1);


set(get(gca,'XLabel'),'String','cumulative number of atoms');
set(get(gca,'YLabel'),'String','cumulative number of excess species');

if file
    saveas(ftmp, [path file]);
    delete(ftmp);
    
else
    plotedit on
    
end
