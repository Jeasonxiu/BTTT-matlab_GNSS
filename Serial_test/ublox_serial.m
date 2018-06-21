function varargout = ublox_serial(varargin)
% UBLOX_SERIAL MATLAB code for ublox_serial.fig
%      UBLOX_SERIAL, by itself, creates a new UBLOX_SERIAL or raises the existing
%      singleton*.
%
%      H = UBLOX_SERIAL returns the handle to a new UBLOX_SERIAL or the handle to
%      the existing singleton*.
%
%      UBLOX_SERIAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UBLOX_SERIAL.M with the given input arguments.
%
%      UBLOX_SERIAL('Property','Value',...) creates a new UBLOX_SERIAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ublox_serial_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ublox_serial_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ublox_serial

% Last Modified by GUIDE v2.5 12-Jul-2017 11:02:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ublox_serial_OpeningFcn, ...
                   'gui_OutputFcn',  @ublox_serial_OutputFcn, ...
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


% --- Executes just before ublox_serial is made visible.
function ublox_serial_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ublox_serial (see VARARGIN)

% Choose default command line output for ublox_serial
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ublox_serial wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ublox_serial_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function lati_Callback(hObject, eventdata, handles)
% hObject    handle to lati (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lati as text
%        str2double(get(hObject,'String')) returns contents of lati as a double


% --- Executes during object creation, after setting all properties.
function lati_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lati (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function longi_Callback(hObject, eventdata, handles)
% hObject    handle to longi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of longi as text
%        str2double(get(hObject,'String')) returns contents of longi as a double


% --- Executes during object creation, after setting all properties.
function longi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to longi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_push.
function run_push_Callback(hObject, eventdata, handles)
% hObject    handle to run_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.run_push,'UserData',1);

%% Com port 설정
comport = get(handles.port, 'String');
COM = strcat('com',comport);

%% Reference Point 지정
ref_point = get(handles.ref, 'String')
VRSx = get(handles.vrsx, 'String');clc
VRSy = get(handles.vrsy, 'String');
VRSz = get(handles.vrsz, 'String');
if ~isempty(ref_point)
    %% 구글맵 로드
    axes(handles.axes3)
    axis([126.8766897178555 126.87760580319954 37.47957614373198 37.480091288616045]);
    hold on
    plot(126.876946215514, 37.4797970132665, 'rd')
    text(126.876946215514, 37.4797970132665, 'A Point','Fontsize',15)
    plot(126.876985221076, 37.4797132062492,' bd')
    text(126.876985221076, 37.4797132062492, 'B Point','Fontsize',15)
    plot_google_map;
    if ref_point =='A' | ref_point =='a'
        TruePos = [-3041235.578 4053941.677 3859881.013]  % : JPspace A point
    elseif ref_point =='B' | ref_point =='b'
        TruePos = [-3041241.741 4053944.143 3859873.640]   % : JPspace B point
    elseif ref_point =='C' | ref_point =='c'
        TruePos = [-3041241.741 4053944.143 3859873.640]   % : JPspace C point
    elseif ref_point =='D' | ref_point =='d'
        TruePos = [-3041241.741 4053944.143 3859873.640]   % : JPspace D point
    else
        TruePos = []
    end
    
elseif ~isempty(VRSx) & ~isempty(VRSy) & ~isempty(VRSz)
    TruePos = [str2num(VRSx), str2num(VRSy), str2num(VRSz)]
    forgoogle = xyz2gd(TruePos);
    axes(handles.axes3)
    axis([forgoogle(2)-0.0005 forgoogle(2)+0.0005 forgoogle(1)-0.0005 forgoogle(1)+0.0005]);
    hold on
    plot_google_map
else
    TruePos = []
end

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

%% error plot 시, circle
cir = [0:1:360];
x1=sind(cir); y1=cosd(cir);

%% Serial data 수신
while (get(handles.run_push,'UserData') ~= 0)
    % Communicating with instrument object, obj1.
    data1 = query(obj1, '*IDN?');
    
    %% NMEA 출력시 구글맵 표출 및 에러 표출
    if length(data1) > 60 & data1(1:6) == '$GPGGA'
        GPGGA{epoch,1} = {data1};
%         disp(data1);
        GGA = char(GPGGA{epoch,1})
   
        EPOCH = EPOCH + 1
%         set(handles.edit20,'String',GGA);
% %         guidata(hObject, handles);
%         if length(data1) > 60
            GGA=GGA(1:68)
            gpgga{EPOCH,1} = {GGA};
            assignin('base','GPGGA',gpgga);
%             [hh,mm,ss,la,lo,qi,ht] = NEWreadGGA2(GGA);
            [hh,mm,ss,x,y,z,la,lo,qi,nSats,ht] = NEWreadGGA3(GGA) 
            set(handles.edit20,'String',GGA);
            set(handles.longi,'String',lo);
            set(handles.lati,'String',la);
            set(handles.height,'String',ht);
            set(handles.qi,'String',qi);
            if ~isempty(TruePos)
                axes(handles.axes4)
                hold on; grid on; 
                plot(0.5*x1, 0.5*y, 'r-')
                plot(x1, y1, 'r:')
                plot(3*x1, 3*y1, 'r-')
                plot(5*x1, 5*y1, 'r-')
                text(0.6,0,'0.5m')
                text(1,1,'1m')
                text(3,1,'3m')
                text(5,1,'5m')
                %% 추정된 XYZ와 참값 XYZ의 차이값 결정
                gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2); 
                dXYZ = [x, y, z] - TruePos;
                %% dXYZ를 dNEV로 변환
                dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
                %% 각 성분별 RMS 계산
                dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
                dNE = sqrt(dN.^2 + dE.^2);        %rmsH = myRMS(dNE);
%                 plot(dNE*x1, dNE*y1,'b:')
                dne = num2str(dNE,'%5.2f');
                set(handles.dNE,'String',dne);
                dv = num2str(dV,'%5.2f');
                set(handles.dV,'String',dv);
                %rmsV = myRMS(dV);
                d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
                d3d = num2str(d3,'%5.2f');
                set(handles.d3D,'String',d3d);
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
            result(EPOCH,:) = [dN, dE, dV, dNE, d3, qi, nSats];
            handles.results = result;
            assignin('base','dNEV',result);
            guidata(hObject, handles);
            axes(handles.axes3)
            hold on;
            plot(lo, la,'r*')
%             set(handles.llh, 'Data', [la, lo, ht]);

            drawnow;
%         end
    elseif length(data1) < 60 & data1(1:6) == '$GPGGA'
        set(handles.edit20,'String',['not fix']);
    else
        GPGGA{epoch,1} ={};
%         assignin('base','GPGGA',GPGGA);
    end
    
    epoch = epoch +1;
%     disp(epoch)
    pause(0.1)
%     set(handles.edit20,'String','');
end

% --- Executes on button press in stop_push.
function stop_push_Callback(hObject, eventdata, handles)
% hObject    handle to stop_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edit20,'String','');
set(handles.run_push,'UserData',0);
stop = get(handles.edit20,'String');
set(handles.edit20,'String',[stop 'STOP pressed.']);
delete(instrfindall)


function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function height_Callback(hObject, eventdata, handles)
% hObject    handle to height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of height as text
%        str2double(get(hObject,'String')) returns contents of height as a double


% --- Executes during object creation, after setting all properties.
function height_CreateFcn(hObject, eventdata, handles)
% hObject    handle to height (see GCBO)
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


% --- Executes on button press in save_push.
function save_push_Callback(hObject, eventdata, handles)
% hObject    handle to save_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
location = get(handles.loca, 'String');
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


function loca_Callback(hObject, eventdata, handles)
% hObject    handle to loca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loca as text
%        str2double(get(hObject,'String')) returns contents of loca as a double


% --- Executes during object creation, after setting all properties.
function loca_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loca (see GCBO)
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


% --- Executes on button press in post_push.
function post_push_Callback(hObject, eventdata, handles)
% hObject    handle to post_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

RESULT = handles.results;
dNE = rms(RESULT(:,4));
d3 = rms(RESULT(:,5));
tHour = [1:1:length(RESULT(:,1))]';

figure();
hold on;
subplot(4,4,[1,2,5,6])
% lineCross(0,0,'r',1)
plot(RESULT(:,2), RESULT(:,1),'bo'); 
axis([-3 3 -3 3]);
axis square
grid on; 
xlabel({'\Delta E (meters)',...
    ['dNE = ', num2str(dNE)],...
    [' 3D = ', num2str(d3)]}); ; ylabel('\Delta N (meters)')
%% 그래프 우측
subplot(4,4,[3,4])
plot(tHour(:,1), RESULT(:,1), '.r:', tHour(:,1), RESULT(:,2), '.b:'); xlim([1 length(tHour)]); grid on;
legend('\Delta N', '\Delta E')
ylabel('\Delta N,E (meters)');
subplot(4,4,[7,8])
plot(tHour(:,1), RESULT(:,3), '.:'); xlim([1 length(tHour)]); grid on;
ylabel('\Delta V (meters)')
legend('\Delta V (meters)')
subplot(4,4,[11,12])
plot(tHour(:,1), RESULT(:,6),'.b:'); xlim([1 length(tHour)]); grid on;
ylabel('Fix Quality');
subplot(4,4,[15,16])
plot(tHour(:,1), RESULT(:,7), '.:'); xlim([1 length(tHour)]); grid on;
ylabel('Number of Satellite')



function vrsx_Callback(hObject, eventdata, handles)
% hObject    handle to vrsx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vrsx as text
%        str2double(get(hObject,'String')) returns contents of vrsx as a double


% --- Executes during object creation, after setting all properties.
function vrsx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vrsx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vrsy_Callback(hObject, eventdata, handles)
% hObject    handle to vrsy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vrsy as text
%        str2double(get(hObject,'String')) returns contents of vrsy as a double


% --- Executes during object creation, after setting all properties.
function vrsy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vrsy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vrsz_Callback(hObject, eventdata, handles)
% hObject    handle to vrsz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vrsz as text
%        str2double(get(hObject,'String')) returns contents of vrsz as a double


% --- Executes during object creation, after setting all properties.
function vrsz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vrsz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dV_Callback(hObject, eventdata, handles)
% hObject    handle to dV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dV as text
%        str2double(get(hObject,'String')) returns contents of dV as a double


% --- Executes during object creation, after setting all properties.
function dV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function d3D_Callback(hObject, eventdata, handles)
% hObject    handle to d3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d3D as text
%        str2double(get(hObject,'String')) returns contents of d3D as a double


% --- Executes during object creation, after setting all properties.
function d3D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
