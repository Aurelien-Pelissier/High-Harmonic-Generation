function varargout = main(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_OpeningFcn, ...
                   'gui_OutputFcn',  @main_OutputFcn, ...
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


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
clc

%inial values
handles.I0 = 7e13;          %peak intensity
handles.V  = 250;           %gas velocity
handles.Ptot = 500;         %maximum pressure
handles.lp = 200e-6;        %interaction length
handles.profile = 'gauss';  %density profile
handles.t = 0;

handles.q		= 21;       %harmonic order
handles.alpha	= 2e-14;    %phase coefficient
handles.tp		= 130e-15;  %pulse duration
handles.lambda1 = 1050e-9;  %fundamental wavelength
handles.R0		= 19.6e-6;  %beam waist
handles.gas     = 'Kr';     %Kr for Krypton and Xe for Xenon
handles.Te = 3;             %freed electron temperature
handles.f = 60e6;           %pulse frequency
handles.xh = 0;       %helium fraction
handles.P = handles.Ptot*(1-handles.xh); %partial pressure of the gas (removing helium)

handles.znozzle = 0;        %nozzle position
handles.zmax	= 1e-3;     %graph parameter zmax
handles.rmax	= 50e-6;    %graph parameter rmax
handles.nres    = 50;       %graph parameter resolution
handles.ioniz_max = 0.1;    %maximum ionization when plot

%GUI parameters (callback and popmenue)
handles.plotpm = 1;
handles.plotio = 1;
handles.plotdi = 1;
handles.plot_noz = 1;
handles.plot_t = 1;
handles.plot_t = 1;
handles.checkz = 1;
handles.checkt = 1;
handles.sav_noz = 0;
handles.sav_t = 1;

%definings variable to check which button has been modified
handles.iinoz1 = 1;
handles.iinoz2 = 1;
handles.iit1 = 1;
handles.iit2 = 1;
handles.iiq = 1;
handles.iiI0 = 1;
handles.iiV = 1;
handles.iiP = 1;
handles.iilambda = 1;
handles.iigas = 1;
handles.iiR0 = 1;
handles.iitp = 1;
handles.iiTe = 1;
handles.iialpha = 1;
handles.iizmax = 1;
handles.iinres = 1;
handles.iif = 1;
handles.iiPow = 1; 
handles.iixh = 1;
handles.iilp = 1;
handles.iiprof = 1;
handles.iiPtot = 1;


%phase matching%
[handles.Dk,handles.Dp,handles.ethai,handles.ethaf,handles.ethafb,handles.ethaft] = phase_matching(handles.I0,handles.V,handles.P,handles.lp,handles.profile,handles.znozzle,handles.gas,handles.q,handles.f,handles.R0,handles.lambda1,handles.tp,handles.Te,handles.alpha,handles.zmax,handles.rmax,handles.nres,handles.t);
handles.Cl = pi./handles.Dk;
handles.current_data_PM = cell(1);
handles.current_data_PM{1} = handles.Cl;
handles.current_data_PM{2} = handles.Dk;
handles.current_data_PM{3} = handles.Dp;

handles.current_data_Io = cell(1);
handles.current_data_Io{1} = handles.ethai;
handles.current_data_Io{2} = handles.ethaf;
handles.current_data_Io{3} = handles.ethafb;
handles.current_data_Io{4} = handles.ethaft;

axes(handles.axes1)
plot_PM(handles.current_data_PM{handles.plotpm},handles.zmax,handles.rmax,handles.nres,handles.znozzle)


axes(handles.axes2)
plot_Io(handles.current_data_Io{handles.plotio},handles.zmax,handles.rmax,handles.nres,handles.ioniz_max)
colormap(handles.axes2,jet)


%dipole response%
[handles.dipq,handles.dipqt,handles.qW,handles.dipo] = detectdipole(handles.I0,handles.lambda1,handles.znozzle,handles.tp,handles.R0,handles.gas,handles.q,handles.t,handles.zmax,handles.nres,handles.axes3);

axes(handles.axes3)
plot_dip(handles.qW,handles.dipo,handles.q,handles.t, handles.zmax,handles.tp,handles.nres,handles.dipq, handles.dipqt, handles.plotdi,handles.znozzle,handles.lp,handles.profile,handles.R0,handles.lambda1)
 

%final amplitude%
    %harmonic growth
    axes(handles.axes4)
    [handles.z,handles.Iq] = amplitude(handles.dipq,handles.ethaf,handles.Dp,handles.P,handles.lp,handles.profile,handles.q,handles.lambda1,handles.znozzle,handles.zmax,handles.nres,handles.gas);
    plot_ampl(handles.Iq,handles.z)

    %nozzle position dependance
%     [handles.Iqz,handles.Clz] = amp_nozzle(handles.dipq,handles.I0,handles.V,handles.P,handles.lp,handles.profile,handles.gas,handles.q,handles.f,handles.R0,handles.lambda1,handles.tp,handles.Te,handles.alpha,handles.zmax,handles.rmax,handles.nres,handles.t);
    handles.Iqz = zeros(handles.nres+1,1); %too long to open
    handles.Clz = zeros(handles.nres+1,1);
    axes(handles.axes5)
    plot_ampnoz(handles.Iqz,handles.Clz,handles.zmax,handles.nres,handles.znozzle,handles.plot_noz,handles.sav_noz)

    %time dependance
%     [handles.Iqt,handles.Clt] = amp_time(handles.I0,handles.V,handles.P,handles.lp,handles.profile,handles.gas,handles.q,handles.f,handles.R0,handles.lambda1,handles.tp,handles.Te,handles.alpha,handles.zmax,handles.rmax,handles.nres,handles.znozzle,handles.axes3);
    handles.Iqt = zeros(handles.nres+1,1); %too long to open
    handles.Clt = zeros(handles.nres+1,1);       
    axes(handles.axes6)
    plot_ampt(handles.Iqt,handles.Clt,handles.tp,handles.nres,handles.t,handles.plot_t)
    

% Choose default command line output for main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)     
    
%phase matching
[handles.Dk,handles.Dp,handles.ethai,handles.ethaf,handles.ethafb,handles.ethaft] = phase_matching(handles.I0,handles.V,handles.P,handles.lp,handles.profile,handles.znozzle,handles.gas,handles.q,handles.f,handles.R0,handles.lambda1,handles.tp,handles.Te,handles.alpha,handles.zmax,handles.rmax,handles.nres,handles.t);
handles.Cl = pi./handles.Dk;
handles.current_data_PM = cell(1);
handles.current_data_PM{1} = handles.Cl;
handles.current_data_PM{2} = handles.Dk;
handles.current_data_PM{3} = handles.Dp;

handles.current_data_Io = cell(1);
handles.current_data_Io{1} = handles.ethai;
handles.current_data_Io{2} = handles.ethaf;
handles.current_data_Io{3} = handles.ethafb;
handles.current_data_Io{4} = handles.ethaft;

axes(handles.axes1)
plot_PM(handles.current_data_PM{handles.plotpm},handles.zmax,handles.rmax,handles.nres,handles.znozzle)


axes(handles.axes2)
plot_Io(handles.current_data_Io{handles.plotio},handles.zmax,handles.rmax,handles.nres,handles.ioniz_max)
colormap(handles.axes2,jet)

%dipole response
[handles.dipq,handles.dipqt,handles.qW,handles.dipo] = detectdipole(handles.I0,handles.lambda1,handles.znozzle,handles.tp,handles.R0,handles.gas,handles.q,handles.t,handles.zmax,handles.nres,handles.axes3);

axes(handles.axes3)
plot_dip(handles.qW,handles.dipo,handles.q,handles.t, handles.zmax,handles.tp,handles.nres,handles.dipq, handles.dipqt, handles.plotdi, handles.znozzle,handles.lp,handles.profile,handles.R0,handles.lambda1)
    
%final amplitude
axes(handles.axes4)
[handles.z,handles.Iq] = amplitude(handles.dipq,handles.ethaf,handles.Dp,handles.P,handles.lp,handles.profile,handles.q,handles.lambda1,handles.znozzle,handles.zmax,handles.nres,handles.gas);
plot_ampl(handles.Iq,handles.z)


%looking for the checkbox and calculating for all znozzle if checked
handles.checkz = get(handles.checkbox1,'value');
productz = handles.iinoz1*handles.iit1*handles.iiq*handles.iiI0*handles.iiV*handles.iiP*handles.iilambda*handles.iigas*handles.iiR0*handles.iitp*handles.iiTe*handles.iialpha*handles.iizmax*handles.iinres*handles.iif*handles.iiPow*handles.iixh*handles.iilp*handles.iiprof*handles.iiPtot;    
if handles.checkz == 1
    if productz == 1
        if any(handles.Iqz) == 0
            [handles.Iqz,handles.Clz] = amp_nozzle(handles.dipq,handles.I0,handles.V,handles.P,handles.lp,handles.profile,handles.gas,handles.q,handles.f,handles.R0,handles.lambda1,handles.tp,handles.Te,handles.alpha,handles.zmax,handles.rmax,handles.nres,handles.t);
        end
    else
        [handles.Iqz,handles.Clz] = amp_nozzle(handles.dipq,handles.I0,handles.V,handles.P,handles.lp,handles.profile,handles.gas,handles.q,handles.f,handles.R0,handles.lambda1,handles.tp,handles.Te,handles.alpha,handles.zmax,handles.rmax,handles.nres,handles.t);
    end
else
    if productz == 0
        handles.Iqz = zeros(handles.nres+1,1);
        handles.Clz = zeros(handles.nres+1,1);
    end
end  


axes(handles.axes5)
plot_ampnoz(handles.Iqz,handles.Clz,handles.zmax,handles.nres,handles.znozzle,handles.plot_noz,handles.sav_noz)


%looking for the checkbox and calculating for all time position if checked
handles.checkt = get(handles.checkbox3,'value');
productt = handles.iinoz2*handles.iit2*handles.iiq*handles.iiI0*handles.iiV*handles.iiP*handles.iilambda*handles.iigas*handles.iiR0*handles.iitp*handles.iiTe*handles.iialpha*handles.iizmax*handles.iinres*handles.iif*handles.iiPow*handles.iixh*handles.iilp*handles.iiprof*handles.iiPtot;
if handles.checkt == 1
    if productt == 1
        if any(handles.Iqt) == 0
            [handles.Iqt,handles.Clt] = amp_time(handles.I0,handles.V,handles.P,handles.lp,handles.profile,handles.gas,handles.q,handles.f,handles.R0,handles.lambda1,handles.tp,handles.Te,handles.alpha,handles.zmax,handles.rmax,handles.nres,handles.znozzle,handles.axes3);
        end
    else
        [handles.Iqt,handles.Clt] = amp_time(handles.I0,handles.V,handles.P,handles.lp,handles.profile,handles.gas,handles.q,handles.f,handles.R0,handles.lambda1,handles.tp,handles.Te,handles.alpha,handles.zmax,handles.rmax,handles.nres,handles.znozzle,handles.axes3);
    end
else
    if productt == 0
        handles.Iqt = zeros(handles.nres+1,1);
        handles.Clt = zeros(handles.nres+1,1);
    end
end 





    axes(handles.axes6)
    plot_ampt(handles.Iqt,handles.Clt,handles.tp,handles.nres,handles.t,handles.plot_t)
    
%reset all the "check if modified" variable
handles.iinoz1 = 1;
handles.iinoz2 = 1;
handles.iit1 = 1;
handles.iit2 = 1;
handles.iiq = 1;
handles.iiI0 = 1;
handles.iiV = 1;
handles.iiP = 1;
handles.iilambda = 1;
handles.iigas = 1;
handles.iiR0 = 1;
handles.iitp = 1;
handles.iiTe = 1;
handles.iialpha = 1;
handles.iizmax = 1;
handles.iinres = 1;
handles.iif = 1;
handles.iiPow = 1; 
handles.iixh = 1;
handles.iilp = 1;
handles.iiprof = 1;
handles.iiPtot = 1;



guidata(hObject, handles);
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

handles.I0 = get(handles.edit1,'string');
handles.I0 = str2double(handles.I0);
handles.iiI0 = 0; 

k1		= 2*pi/handles.lambda1;
Zr		= k1*handles.R0^2/2;
handles.Pow = handles.I0/ (4*sqrt(log(2))/(handles.R0^2*handles.tp*handles.f*pi^(3/2))) *1e4;
set(handles.edit21, 'String' , sprintf('%.1f',handles.Pow*1e-3))
guidata(hObject, handles);


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



function edit2_Callback(hObject, eventdata, handles)

handles.V = get(handles.edit2,'string');
handles.V = str2double(handles.V);
handles.iiV = 0;

np = floor(2.5*handles.R0*handles.f/handles.V)+1;
set(handles.edit23, 'String' , sprintf('%.0f',np))

handles.P = handles.Ptot*(1-handles.xh);
set(handles.edit25, 'String' , sprintf('%.0f',handles.P))

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)

handles.Ptot = get(handles.edit3,'string');
handles.Ptot = str2double(handles.Ptot);
handles.iiPtot = 0; 
handles.P = handles.Ptot*(1-handles.xh);
set(handles.edit25, 'String' , sprintf('%.0f',handles.P))

%absorption length
sigma = absorb(handles.lambda1,handles.q,handles.gas);
Kb = 1.38064852e-23;
rho = handles.P*100/(Kb*293);
handles.Labs = 2/(sigma*rho);
set(handles.edit31, 'String' , sprintf('%.0f',handles.Labs*1e6))

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2



val = get(hObject,'Value');
str = get(hObject,'String');
axes(handles.axes1)
switch str{val}
    case 'Coherence length [m] (log)'
        handles.plotpm = 1;
    case 'DeltaK [m-1] (log)'
        handles.plotpm = 2;
    case 'DeltaPhi'
        handles.plotpm = 3;
end
axes(handles.axes1)
plot_PM(handles.current_data_PM{handles.plotpm},handles.zmax,handles.rmax,handles.nres,handles.znozzle)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit5_Callback(hObject, eventdata, handles)

handles.q = get(handles.edit5,'string');
handles.q = str2double(handles.q);
handles.iiq = 0; 

%absorption length
sigma = absorb(handles.lambda1,handles.q,handles.gas);
Kb = 1.38064852e-23;
rho = handles.P*100/(Kb*293);
handles.Labs = 2/(sigma*rho);
set(handles.edit31, 'String' , sprintf('%.0f',handles.Labs*1e6))

guidata(hObject, handles);


function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)

handles.R0 = get(handles.edit6,'string');
handles.R0 = str2double(handles.R0)*1e-6;;
handles.iiR0 = 0; 

np = floor(2.5*handles.R0*handles.f/handles.V)+1
set(handles.edit23, 'String' , sprintf('%.0f',np))
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)

handles.lambda1 = get(handles.edit7,'string');
handles.lambda1 = str2double(handles.lambda1)*1e-9;
handles.iilambda = 0; 

%absorption length
sigma = absorb(handles.lambda1,handles.q,handles.gas);
Kb = 1.38064852e-23;
rho = handles.P*100/(Kb*293);
handles.Labs = 2/(sigma*rho);
set(handles.edit31, 'String' , sprintf('%.0f',handles.Labs*1e6))

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)

handles.tp = get(handles.edit8,'string');
handles.tp = str2double(handles.tp)*1e-15;
handles.iitp = 0; 
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)

handles.Te = get(handles.edit9,'string');
handles.Te = str2double(handles.Te);
handles.iiTe = 0; 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)

handles.alpha = get(handles.edit10,'string');
handles.alpha = str2double(handles.alpha);
handles.iialpha = 0; 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit14_Callback(hObject, eventdata, handles)

handles.zmax = get(handles.edit14,'string');
handles.zmax = str2double(handles.zmax)*1e-3;
handles.iizmax = 0; 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
handles.nres = get(handles.edit16,'string');
handles.nres = str2double(handles.nres);
handles.iinres = 0; 
guidata(hObject, handles);


function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4

val = get(hObject,'Value');
str = get(hObject,'String');
switch str{val}
    case 'Ionization after 1 pulse'
        handles.plotio = 1;

    case 'Ionization after last pulse' 
        handles.plotio = 2;
        
    case 'Ionization before last pulse (Xe* + Xe+)' 
        handles.plotio = 3; 
        
    case 'Ionization at time position'
        handles.plotio = 4;
end
axes(handles.axes2)
plot_Io(handles.current_data_Io{handles.plotio},handles.zmax,handles.rmax,handles.nres,handles.ioniz_max)
colormap(handles.axes2,jet)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)

val = get(hObject,'Value');
str = get(hObject,'String');
switch str{val}
    case 'Krypton'
        handles.gas = 'Kr';
    case 'Xenon'
        handles.gas = 'Xe';
    case 'Argon'
        handles.gas = 'Ar';
end
handles.iigas = 0; 

%absorption length
sigma = absorb(handles.lambda1,handles.q,handles.gas);
Kb = 1.38064852e-23;
rho = handles.P*100/(Kb*293);
handles.Labs = 2/(sigma*rho);
set(handles.edit31, 'String' , sprintf('%.0f',handles.Labs*1e6))

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)

val = get(hObject,'Value');
str = get(hObject,'String');
axes(handles.axes3)
switch str{val}
    case 'Sw(q) at nozzle and time position'
        handles.plotdi = 1;
    case 'Swq(z) at time position'
        handles.plotdi = 2;           
    case 'Swq(t) at nozzle position'
        handles.plotdi = 3;            
end
axes(handles.axes3)
plot_dip(handles.qW,handles.dipo,handles.q,handles.t, handles.zmax,handles.tp,handles.nres,handles.dipq, handles.dipqt, handles.plotdi, handles.znozzle,handles.lp,handles.profile,handles.R0,handles.lambda1)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8


% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
handles.znozzle = get(handles.edit18,'string');
handles.znozzle = str2double(handles.znozzle)*1e-3;
handles.iinoz1 = 1; 
handles.iinoz2 = 0;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
axes(handles.axes5)
[handles.Iqz,handles.Clz] = amp_nozzle(handles.dipq,handles.I0,handles.V,handles.P,handles.lp,handles.profile,handles.gas,handles.q,handles.f,handles.R0,handles.lambda1,handles.tp,handles.Te,handles.alpha,handles.zmax,handles.rmax,handles.nres,handles.t);
plot_ampnoz(handles.Iqz,handles.Clz,handles.zmax,handles.nres,handles.znozzle,handles.plot_noz,handles.sav_noz)
  
guidata(hObject, handles);

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1



function edit20_Callback(hObject, eventdata, handles)
handles.f = get(handles.edit20,'string');
handles.f = str2double(handles.f)*1e6;
handles.iif = 0; 

np = floor(2.5*handles.R0*handles.f/handles.V)+1;
set(handles.edit23, 'String' , sprintf('%.0f',np))
guidata(hObject, handles);


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
        


function edit21_Callback(hObject, eventdata, handles)
handles.Pow = get(handles.edit21,'string');
handles.Pow = str2double(handles.Pow)*1e3;
handles.iiPow = 0;

k1		= 2*pi/handles.lambda1;
Zr		= k1*handles.R0^2/2;
handles.I0 = handles.Pow *(4*sqrt(log(2))/(handles.R0^2*handles.tp*handles.f*pi^(3/2))) * 1e-4;
set (handles.edit1, 'String' , sprintf('%.1e',handles.I0))
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5



function edit22_Callback(hObject, eventdata, handles)
handles.xh = get(handles.edit22,'string');
handles.xh = str2double(handles.xh)/100;
handles.iixh = 0;

np = floor(2.5*handles.R0*handles.f/handles.V)+1;
set(handles.edit23, 'String' , sprintf('%.0f',np))

handles.P = handles.Ptot*(1-handles.xh);
set(handles.edit25, 'String' , sprintf('%.0f',handles.P))
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6



function edit24_Callback(hObject, eventdata, handles)
handles.lp = get(handles.edit24,'string');
handles.lp = str2double(handles.lp)*1e-6;
handles.iilp = 0;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton8.
function radiobutton8_Callback(hObject, eventdata, handles)
    handles.profile = 'gauss';
    handles.iiprof = 0;
    guidata(hObject, handles);


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
    handles.profile = 'squar';
    handles.iiprof = 0;
    guidata(hObject, handles);
    
    
% --- Executes on button press in radiobutton13.
function radiobutton13_Callback(hObject, eventdata, handles)
    handles.profile = 'squ_g';
    handles.iiprof = 0;
    guidata(hObject, handles);



% --- Executes when selected object is changed in uibuttongroup3.
function uibuttongroup3_SelectionChangedFcn(hObject, eventdata, handles)
    get(eventdata.NewValue, 'Tag');




function edit25_Callback(hObject, eventdata, handles)
handles.P = get(handles.edit25,'string');
handles.P = str2double(handles.P);
handles.iiP = 0;
handles.Ptot = handles.P/(1-handles.xh);
set(handles.edit3, 'String' , sprintf('%.0f',handles.Ptot))

%absorption length
sigma = absorb(handles.lambda1,handles.q,handles.gas);
Kb = 1.38064852e-23;
rho = handles.P*100/(Kb*293);
handles.Labs = 2/(sigma*rho);
set(handles.edit31, 'String' , sprintf('%.0f',handles.Labs*1e6))

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
axes(handles.axes6)
[handles.Iqt,handles.Clt] = amp_time(handles.I0,handles.V,handles.P,handles.lp,handles.profile,handles.gas,handles.q,handles.f,handles.R0,handles.lambda1,handles.tp,handles.Te,handles.alpha,handles.zmax,handles.rmax,handles.nres,handles.znozzle,handles.axes3);
plot_ampt(handles.Iqt,handles.Clt,handles.tp,handles.nres,handles.t,handles.plot_t)

guidata(hObject, handles);


% --- Executes on selection change in popupmenu9.
function popupmenu9_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
str = get(hObject,'String');
switch str{val}
    case 'Harmonic amplitude at nozzle position (t)'
        handles.plot_t = 1;
    case 'Coherence length at nozzle position (t)'
        handles.plot_t = 2;        
end
axes(handles.axes6)
plot_ampt(handles.Iqt,handles.Clt,handles.tp,handles.nres,handles.t,handles.plot_t)

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
str = get(hObject,'String');
switch str{val}
    case 'Harmonic amplitude at time position (znozzle)'
        handles.plot_noz = 1;
    case 'Coherence length at time position (znozzle)'
        handles.plot_noz = 2;           
end
axes(handles.axes5)
plot_ampnoz(handles.Iqz,handles.Clz,handles.zmax,handles.nres,handles.znozzle,handles.plot_noz,handles.sav_noz)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit27_Callback(hObject, eventdata, handles)
handles.t = get(handles.edit27,'string');
handles.t = str2double(handles.t)*1e-15;
handles.iit1 = 0;
handles.iit2 = 1;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit28_Callback(hObject, eventdata, handles)
handles.ioniz_max = get(handles.edit28,'string');
handles.ioniz_max = str2double(handles.ioniz_max)*1e-2;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double


% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit30 as text
%        str2double(get(hObject,'String')) returns contents of edit30 as a double


% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit31_Callback(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit31 as text
%        str2double(get(hObject,'String')) returns contents of edit31 as a double


% --- Executes during object creation, after setting all properties.
function edit31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function plot_PM(data,zmax,rmax,nres,znozzle)
    zmin	= -zmax*(1-0.01); %to avoid z=0 which cause some divergence problem
    rmin	= -rmax;
    z		= [zmin : (zmax-zmin)/nres : zmax];
    r		= [rmin : (rmax-rmin)/nres : rmax];

    
    if max(max(data)) <= 2*pi && min(min(data)) >= 0
        imagesc(z*1e3,r*1e6,data)
        caxis([0,2*pi])
        colormap(parula)
        
    elseif max(max(data)) >= 1e4 
        imagesc(z*1e3,r*1e6,log10(abs(data)))
        colormap(jet)
        caxis([0,6])
        
    else 
        imagesc(z*1e3,r*1e6,log10(abs(data)))
        colormap(jet)
        colormap(flipud(colormap))
        caxis([-6,0])

    end
    
    hold on
    plot([znozzle*1e3,znozzle*1e3],[rmin*1e6,rmax*1e6],'--','linewidth',2,'Color',[0.4,0.4,0.4])
    hold off
    
    set(gca,'XTick',[]);
    ylabel('r[um]')
    colorbar
    
    
function plot_Io(data,zmax,rmax,nres,ioniz_max)
    zmin	= -zmax*(1-0.01); %to avoid z=0 which cause some divergence problem
    rmin	= -rmax;
    z		= [zmin : (zmax-zmin)/nres : zmax];
    r		= [rmin : (rmax-rmin)/nres : rmax];
    imagesc(z.*1e3,r.*1e6,data)
    xlabel('z[mm]'); ylabel('r[um]')
    caxis([0,ioniz_max])
    colorbar
    
    
 function plot_dip(qW,dipo,q,t, zmax,tp,nres,dipqz,dipqt,plotdi,znozzle,lp,profile,R0,lambda1)
    
    if plotdi == 1;
        semilogy(qW,dipo,'linewidth', 1.1,'color','blue');
        hold on
        semilogy([q,q],[1e-20,1e5],'--','linewidth',1,'Color',[0.4,0.4,0.4])
        xlabel('Harmonic order q'); ylabel('S(w) [arb. unit]')
        xlim([1 39])
        ylim([1e-15 1e3])
        grid on
        hold off
    
    elseif plotdi == 2;    
        zmin	= -zmax;
        z		= [zmin : (zmax-zmin)/nres : zmax];
        for i = 1:length(z)
            P(i) = max(dipqz)*Press(z(i),0,1,lp,profile,znozzle);
        end
        plot(z*1e3,dipqz,'linewidth', 2,'color','red');
        hold on
        plot(z*1e3,P,'--','linewidth', 1.5,'Color',[0.4,0.4,0.4])
%         plot([znozzle*1e3,znozzle*1e3],[1.1*min(dipqz),1.1*max(dipqz)],'--','linewidth',1,'Color',[0.4,0.4,0.4])
        hold off
        xlabel('z[mm]'); ylabel('S(w) [arb. unit]')
        ylim([0,1.1*max(dipqz)])
        legend('Dipole response','Density profile')
        
     elseif plotdi == 3; 
        tmesh = [-tp : (2*tp)/nres : tp];
        for i = 1:length(tmesh)
            Ig(i) = max(dipqt)*Igauss(1,tp,R0,lambda1,znozzle,0,tmesh(i));
        end
        
        plot(tmesh*1e15,dipqt,'linewidth', 2,'color','red');
        hold on
        plot([t*1e15,t*1e15],[1.1*min(dipqt),1.1*max(dipqt)],'--','linewidth',1,'Color',[0.4,0.4,0.4])
        plot(tmesh*1e15,Ig,'--','linewidth',1,'Color',[0.6,0.6,0.6])
        hold off
        xlabel('t[fs]'); ylabel('S(w) [arb. unit]')
        ylim([0,1.1*max(dipqt)]); xlim([-tp*1e15,tp*1e15])
    end
    
    
function plot_ampl(E,z)   
    plot(z*1e3,abs(E),'linewidth',2)
    xlabel('z [mm]'), ylabel('Amplitude [arb. unit]')
    if z(end) ~= -z(1)
        xlim([z(1)*1e3,-z(1)*1e3])
    end
    


function [Iqz,Clz] = amp_nozzle(dipqz,I0,V,P,lp,profile,gas,q,f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t)        
      zmin = -zmax;
      znozzle = [zmin : (zmax-zmin)/nres : zmax];
      for  i = 1:length(znozzle)
        [Dk,Dp,~,~,ethaft] = phase_matching(I0,V,P,lp,profile,znozzle(i),gas,q,f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t);
        Clz(i) = pi/Dk(floor(nres/2)+1,i);
        [~,Iq] = amplitude(dipqz,ethaft,Dp,P,lp,profile,q,lambda1,znozzle(i),zmax,nres,gas);
        Iqz(i) = Iq(end);
      end
      
function [Iqt,Clt] = amp_time(I0,V,P,lp,profile,gas,q,f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,znozzle,axes3)        
      tmesh = [-tp : 2*tp/nres : tp];
      for  i = 1:length(tmesh)
        dipqz = detectdipole(I0,lambda1,znozzle,tp,R0,gas,q,tmesh(i),zmax,nres,axes3);
        [Dk,Dp,~,~,ethaft] = phase_matching(I0,V,P,lp,profile,znozzle,gas,q,f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,tmesh(i));
        Clt(i) = pi/Dk(floor(nres/2)+1,i);
        [~,Iq] = amplitude(dipqz,ethaft,Dp,P,lp,profile,q,lambda1,znozzle,zmax,nres,gas);
        Iqt(i) = Iq(end);
      end
      
function plot_ampnoz(Iqz,Clz,zmax,nres,znozzle,plot_noz,sav_noz)
    zmin = -zmax;
    z = [zmin : (zmax-zmin)/nres : zmax];    
    if plot_noz == 1;
        plot(z*1e3,abs(Iqz),'linewidth',2,'color','green')
        hold on
        plot([znozzle*1e3,znozzle*1e3],[0.9*min(abs(Iqz)),1.1*max(abs(Iqz))],'--','linewidth',1,'Color',[0.4,0.4,0.4])
        hold off
        if any(Iqz) ~= 0
            ylim([0.9*min(abs(Iqz)),1.1*max(abs(Iqz))])
        end
        ylabel('Amplitude [arb. unit]')
        
    elseif plot_noz == 2;
        semilogy(z*1e3,abs(Clz),'linewidth',2,'color','cyan')
        hold on
        plot([znozzle*1e3,znozzle*1e3],[0.9*min(abs(Clz)),1.1*max(abs(Clz))],'--','linewidth',1,'Color',[0.4,0.4,0.4])
        hold off
        if any(Clz) ~= 0
            ylim([0.9*min(abs(Clz)),1.1*max(abs(Clz))])
        end
        ylabel('Coherence length [m]')
    end
    xlabel('nozzle position [mm]'); 
    
    
function plot_ampt(Iqt,Clt,tp,nres,t,plot_t)
    tmesh = [-tp : (2*tp)/nres : tp];
    if plot_t == 1;
        for i = 1:length(tmesh)
            I(i) = max(abs(Iqt))*Igauss(1,tp,1,1,0,0,tmesh(i));
        end
        plot(tmesh*1e15,abs(Iqt),'linewidth',2,'color','green')
        hold on
        plot(tmesh*1e15,I,'--','linewidth', 1.5,'Color',[0.4,0.4,0.4])
        plot([t*1e15,t*1e15],[0.9*min(abs(Iqt)),1.1*max(abs(Iqt))],'--','linewidth',1,'Color',[0.4,0.4,0.4])
        hold off
        if any(Iqt) ~= 0
            ylim([0.9*min(abs(Iqt)),1.1*max(abs(Iqt))])
        end
        ylabel('Amplitude [arb. unit]')
%         legend('Harmonic signal','Laser pulse')
        
    elseif plot_t == 2;
        for i = 1:length(tmesh)
            I(i) = max(abs(Clt))*Igauss(1,tp,1,1,0,0,tmesh(i));
        end
        semilogy(tmesh*1e15,abs(Clt),'linewidth',2,'color','cyan')
        hold on
        plot(tmesh*1e15,I,'--','linewidth', 1.5,'Color',[0.4,0.4,0.4])
        plot([t*1e15,t*1e15],[0.9*min(abs(Clt)),1.1*max(abs(Clt))],'--','linewidth',1,'Color',[0.4,0.4,0.4])
        hold off
        if any(Clt) ~= 0
            ylim([0.9*min(abs(Clt)),1.1*max(abs(Clt))])
        end
        ylabel('Coherence length [m]')
%         legend('Coherence length','Laser pulse')
    end
    xlabel('time [fs]'); 
    xlim([-tp*1e15, tp*1e15])

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
    %create a time movie
        tmesh = [-handles.tp : (2*handles.tp)/handles.nres : handles.tp];
        zmin	= -handles.zmax*(1-0.01); %to avoid z=0 which cause some divergence problem
        rmin	= -handles.rmax;
        z		= [zmin : (handles.zmax-zmin)/handles.nres : handles.zmax];
        r		= [rmin : (handles.rmax-rmin)/handles.nres : handles.rmax];
        for i = 1:length(tmesh)
            I(i) = max(abs(handles.Clt))*Igauss(1,handles.tp,1,1,0,0,tmesh(i));
        end
        
        for  i = 1:length(tmesh)
            [Dk,~,~,~,~,ethaft] = phase_matching(handles.I0,handles.V,handles.P,handles.lp,handles.profile,handles.znozzle,handles.gas,handles.q,handles.f,handles.R0,handles.lambda1,handles.tp,handles.Te,handles.alpha,handles.zmax,handles.rmax,handles.nres,tmesh(i));
            Cl = pi./Dk;
        
            fig = figure('visible','off');
            h1 = subplot(4,6,[1,2,3,7,8,9]);
                imagesc(z*1e3,r*1e6,log10(abs(Cl)))
                colormap(jet)
                colormap(flipud(colormap))
                caxis([-6,0])
                hold on
                plot([handles.znozzle*1e3,handles.znozzle*1e3],[rmin*1e6,handles.rmax*1e6],'--','linewidth',2,'Color',[0.4,0.4,0.4])
                hold off
                xlabel('z[mm]');
                ylabel('r[um]')
                
            h2 = subplot(4,6,[4,5,6,10,11,12]);
                imagesc(z.*1e3,r.*1e6,ethaft)
                xlabel('z[mm]');
                caxis([0,handles.ioniz_max])
                colormap(h2,jet)
                set(gca,'YtickLabel',[]);
                
                txt1 = sprintf('%.0ffs',tmesh(i)*1e15);
                text(handles.zmax*1e3/2,rmin*1e6*5/6,txt1,'FontSize',12,'color','white')
                        
            h3 = subplot(4,6,[13,14,15,16,17,18]);
                semilogy(tmesh*1e15,abs(handles.Clt),'linewidth',2,'color','blue')
                hold on
                semilogy(tmesh*1e15,I,'--','linewidth', 1.5,'Color',[0.4,0.4,0.4])
                semilogy([tmesh(i)*1e15,tmesh(i)*1e15],[1e-15,2*max(abs(handles.Clt))],'--','linewidth',1,'Color',[0.4,0.4,0.4])
                hold off
                ylabel('LCoh[m]')
                set(gca,'XtickLabel',[]);
                xlim([-handles.tp*1e15, handles.tp*1e15])
                if any(handles.Clt) ~= 0
                    ylim([0.9*min(abs(handles.Clt)),1.1*max(abs(handles.Clt))])
                end
                
            h4 = subplot(4,6,[19,20,21,22,23,24]);
                plot(tmesh*1e15,abs(handles.Iqt),'linewidth',2,'color','red')
                hold on
                plot(tmesh*1e15,I,'--','linewidth', 1.5,'Color',[0.4,0.4,0.4])
                plot([tmesh(i)*1e15,tmesh(i)*1e15],[0,2*max(abs(handles.Iqt))],'--','linewidth',1,'Color',[0.4,0.4,0.4])
                hold off
                ylabel('Ha.amp[arb.u]')
                xlabel('time [fs]'); 
                xlim([-handles.tp*1e15, handles.tp*1e15])
                if any(handles.Iqt) ~= 0
                    ylim([0,1.1*max(abs(handles.Iqt))])
                end
        
            mov(i) = getframe(fig);
        end
        movie2avi(mov, 'phase_matching_t.avi', 'compression', 'None');
    


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)    
    %create a nozzle movie
    znozzle = [-handles.zmax : (2*handles.zmax)/handles.nres : handles.zmax];
        zmin	= -handles.zmax*(1-0.01); %to avoid z=0 which cause some divergence problem
        rmin	= -handles.rmax;
        z		= [zmin : (handles.zmax-zmin)/handles.nres : handles.zmax];
        r		= [rmin : (handles.rmax-rmin)/handles.nres : handles.rmax];
        for i = 1:length(znozzle)
            I(i) = max(abs(handles.Clz))*Igauss(1,handles.tp,1,1,0,0,handles.t);
        end
        
        for  i = 1:length(znozzle)
            [Dk,~,~,~,~,ethaft] = phase_matching(handles.I0,handles.V,handles.P,handles.lp,handles.profile,znozzle(i),handles.gas,handles.q,handles.f,handles.R0,handles.lambda1,handles.tp,handles.Te,handles.alpha,handles.zmax,handles.rmax,handles.nres,handles.t);
            Cl = pi./Dk;
        
            fig = figure('visible','off');
            h1 = subplot(4,6,[1,2,3,7,8,9]);
                imagesc(z*1e3,r*1e6,log10(abs(Cl)))
                colormap(jet)
                colormap(flipud(colormap))
                caxis([-6,0])
                hold on
                plot([znozzle(i)*1e3,znozzle(i)*1e3],[rmin*1e6,handles.rmax*1e6],'--','linewidth',2,'Color',[0.4,0.4,0.4])
                hold off
                xlabel('z[mm]');
                ylabel('r[um]')
                
            h2 = subplot(4,6,[4,5,6,10,11,12]);
                imagesc(z.*1e3,r.*1e6,ethaft)
                xlabel('z[mm]');
                caxis([0,handles.ioniz_max])
                hold on
                plot([znozzle(i)*1e3,znozzle(i)*1e3],[rmin*1e6,handles.rmax*1e6],'--','linewidth',2,'Color',[0.4,0.4,0.4])
                hold off
                set(gca,'YtickLabel',[]);
                colormap(h2,jet)
                
                txt1 = sprintf('%.2fmm',znozzle(i)*1e3);
                text(handles.zmax*1e3/3,rmin*1e6*5/6,txt1,'FontSize',12,'color','white')
                        
            h3 = subplot(4,6,[13,14,15,16,17,18]);
                semilogy(znozzle*1e3,abs(handles.Clz),'linewidth',2,'color','blue')
                hold on
                semilogy([znozzle(i)*1e3,znozzle(i)*1e3],[1e-15,2*max(abs(handles.Clz))],'--','linewidth',1,'Color',[0.4,0.4,0.4])
                hold off
                ylabel('LCoh[m]') 
                set(gca,'XtickLabel',[]);
                xlim([-handles.zmax*1e3, handles.zmax*1e3])
                if any(handles.Clz) ~= 0
                    ylim([0.9*min(abs(handles.Clz)),1.1*max(abs(handles.Clz))])
                end
                
             h4 = subplot(4,6,[19,20,21,22,23,24]);
                plot(znozzle*1e3,abs(handles.Iqz),'linewidth',2,'color','red')
                hold on
                plot([znozzle(i)*1e3,znozzle(i)*1e3],[1e-4,2*max(abs(handles.Iqz))],'--','linewidth',1,'Color',[0.4,0.4,0.4])
                hold off
                ylabel('Ha.amp[arb.u]')
                xlabel('nozzle position [mm]'); 
                xlim([-handles.zmax*1e3, handles.zmax*1e3])
                if any(handles.Iqz) ~= 0
                    ylim([0,1.1*max(abs(handles.Iqz))])
                end
                
        
            mov(i) = getframe(fig);
        end
        movie2avi(mov, 'phase_matching_znozzle.avi', 'compression', 'None');
