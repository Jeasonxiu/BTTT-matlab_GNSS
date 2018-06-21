function varargout = RTAP_serial_test(varargin)
% RTAP_SERIAL_TEST MATLAB code for RTAP_serial_test.fig
%      RTAP_SERIAL_TEST, by itself, creates a new RTAP_SERIAL_TEST or raises the existing
%      singleton*.
%
%      H = RTAP_SERIAL_TEST returns the handle to a new RTAP_SERIAL_TEST or the handle to
%      the existing singleton*.
%
%      RTAP_SERIAL_TEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RTAP_SERIAL_TEST.M with the given input arguments.
%
%      RTAP_SERIAL_TEST('Property','Value',...) creates a new RTAP_SERIAL_TEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RTAP_serial_test_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RTAP_serial_test_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RTAP_serial_test

% Last Modified by GUIDE v2.5 26-Oct-2016 08:19:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @RTAP_serial_test_OpeningFcn, ...
    'gui_OutputFcn',  @RTAP_serial_test_OutputFcn, ...
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


% --- Executes just before RTAP_serial_test is made visible.
function RTAP_serial_test_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RTAP_serial_test (see VARARGIN)

% Choose default command line output for RTAP_serial_test
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RTAP_serial_test wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RTAP_serial_test_OutputFcn(hObject, eventdata, handles)
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

%% Com port 설정
comport = get(handles.port, 'String');
COM = strcat('com',comport);

%% AP SSID, PW 설정
ssid = get(handles.ssid,'String');
pw = get(handles.pw,'String');
APinfo = strcat('#AP,',ssid,',',pw);

%% Reference Point, RTAP mode 지정
ref_point = get(handles.ref, 'String');
if ref_point =='A' | ref_point =='a'
    TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
elseif ref_point =='B' | ref_point =='b'
    TruePos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
else
    TruePos = [];
end

MODE = get(handles.Mode, 'String');
set_mode = strcat('#mode, ', MODE);

%% serial port 연결 설정
obj1=serial(COM,'BaudRate',115200,'DataBits',8,'StopBits' ,1);

set(obj1,'Parity','none');
warning('off','MATLAB:serial:fscanf:unsuccessfulRead');
% Create the serial port object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = serial(COM);
else
    fclose(obj1);
    obj1 = obj1(1)
end

% Connect to instrument object, obj1.
fopen(obj1);
epoch = 1;
EPOCH = 0;

%% 구글맵 로드
axes(handles.axes3)
axis([126.8766897178555 126.87760580319954 37.47957614373198 37.480091288616045]);
hold on
plot(126.876946215514, 37.4797970132665, 'rd')
text(126.876946215514, 37.4797970132665, 'A Point','Fontsize',15)
plot(126.876985221076, 37.4797132062492,' bd')
text(126.876985221076, 37.4797132062492, 'B Point','Fontsize',15)
plot_google_map;

%% Serial data 수신
while (get(handles.pushbutton1,'UserData') ~= 0)
    % Communicating with instrument object, obj1.
    data1 = query(obj1, '*IDN?');
    %% RTAP 단말기 부팅 완료시 동작 메시지 입력
    if length(data1) >=5 & data1(1:5) == 'Start'
%         fprintf(obj1, '%s\n','#AP,jpspace,43054305')
        fprintf(obj1, '%s\n',APinfo)
        DataToSend = '0xd';
        fprintf(obj1,'%x', DataToSend);
        fprintf(obj1, '%s\n',set_mode)
    end
    
    %% NMEA 출력시 구글맵 표출 및 에러 표출
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
%             [hh,mm,ss,la,lo,qi,ht] = NEWreadGGA2(GGA);
            [hh,mm,ss,x,y,z,la,lo,qi,ht] = NEWreadGGA(GGA) 
            set(handles.edit1,'String',GGA);
            set(handles.lo,'String',lo);
            set(handles.la,'String',la);
            set(handles.ht,'String',ht);
            set(handles.qi,'String',qi);
            if ~isempty(TruePos)
                axes(handles.axes4)
                hold on; grid on; 
                %% 추정된 XYZ와 참값 XYZ의 차이값 결정
                gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2); 
                dXYZ = [x, y, z] - TruePos;
                %% dXYZ를 dNEV로 변환
                dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
                %% 각 성분별 RMS 계산
                dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
                dNE = sqrt(dN.^2 + dE.^2);        %rmsH = myRMS(dNE);
                dne = num2str(dNE,'%5.2f');
                set(handles.dNE,'String',dne);
                %rmsV = myRMS(dV);
                d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
                rXY = max([abs(min(dNEV(:,1))) abs(max(dNEV(:,1))) abs(min(dNEV(:,2))) abs(max(dNEV(:,2)))]);
                rXY = ceil(rXY);
                axis([-rXY rXY -rXY rXY]); grid on; 
                xlabel('\Delta E (meters)'); ylabel('\Delta N (meters)')
                lineCross(0,0,'r',1)
                plot(dE, dN,'ro'); 
                drawnow;
            end
            handles.hms = [hh, mm, ss];
            handles.savegga = gpgga;
            guidata(hObject, handles);
            axes(handles.axes3)
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



function ref_Callback(hObject, eventdata, handles)
% hObject    handle to ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ref as text
%        str2double(get(hObject,'String')) returns contents of ref as a double


% --- Executes during object creation, after setting all properties.
function ref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mode_Callback(hObject, eventdata, handles)
% hObject    handle to Mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mode as text
%        str2double(get(hObject,'String')) returns contents of Mode as a double


% --- Executes during object creation, after setting all properties.
function Mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dNE_Callback(hObject, eventdata, handles)
% hObject    handle to dNE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dNE as text
%        str2double(get(hObject,'String')) returns contents of dNE as a double


% --- Executes during object creation, after setting all properties.
function dNE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dNE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function port_Callback(hObject, eventdata, handles)
% hObject    handle to port (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of port as text
%        str2double(get(hObject,'String')) returns contents of port as a double


% --- Executes during object creation, after setting all properties.
function port_CreateFcn(hObject, eventdata, handles)
% hObject    handle to port (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ssid_Callback(hObject, eventdata, handles)
% hObject    handle to ssid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ssid as text
%        str2double(get(hObject,'String')) returns contents of ssid as a double


% --- Executes during object creation, after setting all properties.
function ssid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ssid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pw_Callback(hObject, eventdata, handles)
% hObject    handle to pw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pw as text
%        str2double(get(hObject,'String')) returns contents of pw as a double


% --- Executes during object creation, after setting all properties.
function pw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




