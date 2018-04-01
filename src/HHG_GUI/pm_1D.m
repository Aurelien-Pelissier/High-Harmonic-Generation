function varargout = pm_1D(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pm_1D_OpeningFcn, ...
                   'gui_OutputFcn',  @pm_1D_OutputFcn, ...
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


% --- Executes just before pm_1D is made visible.
function pm_1D_OpeningFcn(hObject, eventdata, handles, varargin)
clc

handles.I0 = 7e13; 
handles.V  = 300;
handles.Pmax = 2500;
handles.lp = 200e-6;        %interaction length
handles.profile = 'gauss';  %density profile

handles.q		= 21;       %harmonic order
handles.alpha	= 2e-14;    %phase coefficient
handles.tp		= 130e-15;  %pulse duration
handles.lambda1 = 1050e-9;  %fundamental wavelength
handles.R0		= 19.6e-6;  %beam waist
handles.gas     = 'Kr';     %Kr for Krypton and Xe for Xenon
handles.Te = 3;             %freed electron temperature
handles.f = 60e6;           %pulse frequency

handles.znozzle = 0;
handles.zmax	= 0.5e-3;
handles.rmax	= 50e-6;
handles.nres    = 50; %resolution

handles.approx = 1;
handles.plotpm = 1;



[handles.Dk,handles.Dp,handles.ethai,handles.ethaf] = phase_matching_1D_pressure(handles.I0,handles.V,handles.Pmax,handles.lp,handles.profile,handles.gas,handles.q,handles.f,handles.R0,handles.lambda1,handles.tp,handles.Te,handles.alpha,handles.zmax,handles.nres,handles.znozzle,handles.approx);
handles.Cl = pi./handles.Dk;
handles.current_data_PM = cell(1);
handles.current_data_PM{1} = handles.Cl;
handles.current_data_PM{2} = handles.Dk;
handles.current_data_PM{3} = handles.Dp;
handles.current_data_Io{1} = handles.ethai;
handles.current_data_Io{2} = handles.ethaf;

axes(handles.axes1)
plot_PM(handles.current_data_PM{1},handles.zmax,handles.Pmax,handles.znozzle,handles.nres)


% Choose default command line output for pm_1D
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pm_1D wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pm_1D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

[handles.Dk,handles.Dp,handles.ethai,handles.ethaf] = phase_matching_1D_pressure(handles.I0,handles.V,handles.Pmax,handles.lp,handles.profile,handles.gas,handles.q,handles.f,handles.R0,handles.lambda1,handles.tp,handles.Te,handles.alpha,handles.zmax,handles.nres,handles.znozzle,handles.approx);
handles.Cl = pi./handles.Dk;
handles.current_data_PM = cell(1);
handles.current_data_PM{1} = handles.Cl;
handles.current_data_PM{2} = handles.Dk;
handles.current_data_PM{3} = handles.Dp;
handles.current_data_Io{1} = handles.ethai;
handles.current_data_Io{2} = handles.ethaf;

axes(handles.axes1)
plot_PM(handles.current_data_PM{handles.plotpm},handles.zmax,handles.Pmax,handles.znozzle,handles.nres)
guidata(hObject, handles);



% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)

val = get(hObject,'Value');
str = get(hObject,'String');
axes(handles.axes1)
switch str{val}
    case 'Coherence length [m]'
        handles.plotpm = 1;
    case 'DeltaK [m-1]'
        handles.plotpm = 2;
    case 'DeltaPhi'
        handles.plotpm = 3;
end
plot_PM(handles.current_data_PM{handles.plotpm},handles.zmax,handles.Pmax,handles.znozzle,handles.nres)
guidata(hObject, handles);

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



function edit13_Callback(hObject, eventdata, handles)
handles.zmax = get(handles.edit13,'string');
handles.zmax = str2double(handles.zmax)*1e-3;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
handles.nres = get(handles.edit14,'string');
handles.nres = str2double(handles.nres);
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



function edit7_Callback(hObject, eventdata, handles)
handles.q = get(handles.edit7,'string');
handles.q = str2double(handles.q);
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
handles.R0 = get(handles.edit8,'string');
handles.R0 = str2double(handles.R0)*1e-6;;
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
handles.lambda1 = get(handles.edit9,'string');
handles.lambda1 = str2double(handles.lambda1)*1e-9;
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
handles.tp = get(handles.edit10,'string');
handles.tp = str2double(handles.tp)*1e-15;
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
handles.Te = get(handles.edit11,'string');
handles.Te = str2double(handles.Te);
guidata(hObject, handles);


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



function edit12_Callback(hObject, eventdata, handles)
handles.alpha = get(handles.edit12,'string');
handles.alpha = str2double(handles.alpha);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
handles.I0 = get(handles.edit4,'string');
handles.I0 = str2double(handles.I0);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
handles.V = get(handles.edit5,'string');
handles.V = str2double(handles.V);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
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


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
str = get(hObject,'String');
axes(handles.axes1)
switch str{val}
    case 'Approched Calculation (much faster)'
        handles.approx = 1;
    case 'Complete Calculation (recalculate Ionization for each pressure)'
        handles.approx = 0;
end

[handles.Dk,handles.Dp,handles.ethai,handles.ethaf] = phase_matching_1D_pressure(handles.I0,handles.V,handles.Pmax,handles.lp,handles.profile,handles.gas,handles.q,handles.f,handles.R0,handles.lambda1,handles.tp,handles.Te,handles.alpha,handles.zmax,handles.nres,handles.znozzle,handles.approx);
handles.Cl = pi./handles.Dk;
handles.current_data_PM = cell(1);
handles.current_data_PM{1} = handles.Cl;
handles.current_data_PM{2} = handles.Dk;
handles.current_data_PM{3} = handles.Dp;
handles.current_data_Io{1} = handles.ethai;
handles.current_data_Io{2} = handles.ethaf; 
        
plot_PM(handles.current_data_PM{handles.plotpm},handles.zmax,handles.Pmax,handles.znozzle,handles.nres)
        
guidata(hObject, handles);
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
    


function edit18_Callback(hObject, eventdata, handles)
handles.Pmax = get(handles.edit18,'string');
handles.Pmax = str2double(handles.Pmax);
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



function edit20_Callback(hObject, eventdata, handles)
handles.znozzle = get(handles.edit20,'string');
handles.znozzle = str2double(handles.znozzle)*1e-3;
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
handles.f = get(handles.edit21,'string');
handles.f = str2double(handles.f);
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


function plot_PM(data,zmax,Pmax,znozzle,nres)
    zmin	= -zmax*(1-0.01); %to avoid z=0 which cause some divergence problem
    Pmin	= 0;
    z		= [zmin : (zmax-zmin)/nres : zmax];
    P		= [Pmin : (Pmax-Pmin)/nres : Pmax];

    if max(max(data)) >= 2*pi*0.99 && max(max(data)) <= 2*pi*1.01
        imagesc(z*1e3,P,data)
        caxis([0,2*pi])
        colormap(parula)
 
        
    elseif max(max(data)) >= 1e2 
        imagesc(z*1e3,P,log10(abs(data)))
        colormap(jet)
        caxis([0,6]) 
        
    else 
        imagesc(z*1e3,P,log10(abs(data)))
        colormap(jet)
        colormap(flipud(colormap))
        caxis([-6,0]) 

    end
    
    axis xy
    hold on
    plot([znozzle*1e3,znozzle*1e3],[0,Pmax],'--','linewidth',2,'Color',[0.4,0.4,0.4])
    hold off
    
    xlabel('z[mm]'); ylabel('P[mbar]')
    colorbar



function edit22_Callback(hObject, eventdata, handles)
handles.lp = get(handles.edit22,'string');
handles.lp = str2double(handles.lp)*1e-6;
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
