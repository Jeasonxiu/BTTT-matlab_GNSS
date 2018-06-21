function varargout = RTAP_serial_test(varargin)
% GUIGPS MATLAB code for guigps.fig
%      GUIGPS, by itself, creates a new GUIGPS or raises the existing
%      singleton*.
%
%      H = GUIGPS returns the handle to a new GUIGPS or the handle to
%      the existing singleton*.
%
%      GUIGPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIGPS.M with the given input arguments.
%
%      GUIGPS('Property','Value',...) creates a new GUIGPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guigps_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guigps_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guigps

% Last Modified by GUIDE v2.5 25-Oct-2016 10:46:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @guigps_OpeningFcn, ...
    'gui_OutputFcn',  @guigps_OutputFcn, ...
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


% --- Executes just before guigps is made visible.
function guigps_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guigps (see VARARGIN)

% Choose default command line output for guigps
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guigps wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guigps_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton1,'UserData',1);

%% serial port 연결 설정
obj1=serial('COM11','BaudRate',115200,'DataBits',8,'StopBits' ,1);
set(obj1,'Parity','none');
warning('off','MATLAB:serial:fscanf:unsuccessfulRead');
% Create the serial port object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = serial('COM11');
else
    fclose(obj1);
    obj1 = obj1(1)
end

% Connect to instrument object, obj1.
fopen(obj1);
epoch = 1;
EPOCH = 0;
axis([126.87676935370362 126.87736279836793 37.47976946251103 37.48010317548329]);
plot_google_map;
% hold on
while (get(handles.pushbutton1,'UserData') ~= 0)
    % Communicating with instrument object, obj1.
    data1 = query(obj1, '*IDN?');
    if length(data1) >=5 & data1(1:5) == 'Start'
        fprintf(obj1, '%s\n','#AP,jpspace,43054305')
        DataToSend = '0xd';
        fprintf(obj1,'%x', DataToSend);
    end
    
    if length(data1) > 60 & data1(1:6) == '$GPGGA'
        GPGGA{epoch,1} = {data1};
%         disp(data1);
        GGA = char(GPGGA{epoch,1})
   
        EPOCH = EPOCH + 1
%         set(handles.edit1,'String',GGA);
% %         guidata(hObject, handles);
%         if length(data1) > 60
            GGA=GGA(1:68)
            gpgga{EPOCH,1} = {GGA};
            assignin('base','GPGGA',gpgga);
            [hh,mm,ss,la,lo,qi,ht] = NEWreadGGA2(GGA)
            set(handles.edit1,'String',GGA);
            set(handles.lo,'String',lo);
            set(handles.la,'String',la);
            set(handles.ht,'String',ht);
            set(handles.qi,'String',qi);
            handles.hms = [hh, mm, ss];
            handles.savegga = gpgga;
            guidata(hObject, handles);
            hold on;
            plot(lo, la,'r*')
%             set(handles.llh, 'Data', [la, lo, ht]);

            drawnow;
%         end
    else
        GPGGA{epoch,1} ={};
%         assignin('base','GPGGA',GPGGA);
    end
    
    epoch = epoch +1;
%     disp(epoch)
    pause(0.1)
%     set(handles.edit1,'String','');
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(handles.edit1,'String','');
set(handles.pushbutton1,'UserData',0);
stop = get(handles.edit1,'String');
set(handles.edit1,'String',[stop 'STOP pressed.']);
delete(instrfindall)



function la_Callback(hObject, eventdata, handles)
% hObject    handle to la (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of la as text
%        str2double(get(hObject,'String')) returns contents of la as a double


% --- Executes during object creation, after setting all properties.
function la_CreateFcn(hObject, eventdata, handles)
% hObject    handle to la (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lo_Callback(hObject, eventdata, handles)
% hObject    handle to lo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lo as text
%        str2double(get(hObject,'String')) returns contents of lo as a double


% --- Executes during object creation, after setting all properties.
function lo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ht_Callback(hObject, eventdata, handles)
% hObject    handle to ht (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ht as text
%        str2double(get(hObject,'String')) returns contents of ht as a double


% --- Executes during object creation, after setting all properties.
function ht_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ht (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function qi_Callback(hObject, eventdata, handles)
% hObject    handle to qi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of qi as text
%        str2double(get(hObject,'String')) returns contents of qi as a double


% --- Executes during object creation, after setting all properties.
function qi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to qi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
location = get(handles.Loca, 'String');
DOY = get(handles.doy, 'String');
HMS = handles.hms;
GPGGA = handles.savegga;
hour = HMS(1);
switch hour
    case 0
        HOUR = 'a';
    case 1
        HOUR = 'b';
    case 2
        HOUR = 'c';
    case 3
        HOUR = 'd';
    case 4
        HOUR = 'e';
    case 5
        HOUR = 'f';
    case 6
        HOUR = 'g';
    case 7
        HOUR = 'h';
    case 8
        HOUR = 'i';
    case 9
        HOUR = 'j';
    case 10
        HOUR = 'k';
    case 11
        HOUR = 'l';
    case 12
        HOUR = 'm';
    case 13
        HOUR = 'n';
    case 14
        HOUR = 'o';
    case 15
        HOUR = 'p';
    case 16
        HOUR = 'q';
    case 17
        HOUR = 'r';
    case 18
        HOUR = 's';
    case 19
        HOUR = 't';
    case 20
        HOUR = 'u';
    case 21
        HOUR = 'v';
    case 22
        HOUR = 'w';
    case 23
        HOUR = 'x';
end
Filename = strcat(location, DOY, HOUR,'.txt')
file_output = fopen(Filename, 'w');
% gpgga = evalin('base','GPGGA');
for i = 1:length(GPGGA(:,1))
    line = char(GPGGA{i,1});
    fprintf(file_output, '%s\r\n', line);
end



function Loca_Callback(hObject, eventdata, handles)
% hObject    handle to Loca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Loca as text
%        str2double(get(hObject,'String')) returns contents of Loca as a double


% --- Executes during object creation, after setting all properties.
function Loca_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Loca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function doy_Callback(hObject, eventdata, handles)
% hObject    handle to doy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of doy as text
%        str2double(get(hObject,'String')) returns contents of doy as a double


% --- Executes during object creation, after setting all properties.
function doy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to doy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
