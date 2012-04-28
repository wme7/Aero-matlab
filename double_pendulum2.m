function varargout = double_pendulum2(varargin)

% the function simulate double pendulum, using 4th-order Runge-Kutta algorithm for the diffrential equations.
% the diffrential equation are very similar to eq. (12) and (13) here:
%                           http://home1.fvcc.edu/~dhicketh/DiffEqns/spring09projects/LauraStickel/Double%20Pendulum.pdf
% where omega=d(theta)/dt
% Moshe Lindner & Orit Peleg, August 2010 (C)




% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @double_pendulum_OpeningFcn, ...
                   'gui_OutputFcn',  @double_pendulum_OutputFcn, ...
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


% --- Executes just before double_pendulum is made visible.
function double_pendulum_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to double_pendulum (see VARARGIN)

% Choose default command line output for double_pendulum
handles.output = hObject;
handles.len1=.5;
handles.len2=.5;
handles.mass1=.5;
handles.mass2=.5;
handles.t1=0;
handles.t2=0;
handles.w1=0;
handles.w2=0;
handles.G=9.8;
pend_plot(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes double_pendulum wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = double_pendulum_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in start_stop.
function start_stop_Callback(hObject, eventdata, handles)
% hObject    handle to start_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'string','Stop');
set(handles.omega1,'enable','off');
set(handles.theta1,'enable','off');
set(handles.omega2,'enable','off');
set(handles.theta2,'enable','off');
set(handles.length_rat,'enable','off');
set(handles.mass_rat,'enable','off');
set(handles.reset,'enable','off');
set(handles.gravity,'enable','off');
i=1;
while get(hObject,'value')==1
    if (get(handles.traj,'value'))==1
        x_1=handles.len1*sin(handles.t1);
        y_1=-handles.len1*cos(handles.t1) ;
        x_2=handles.len2*sin(handles.t2) + x_1;
        y_2=-handles.len2*cos(handles.t2) + y_1;
        handles.trajectory(1,i)=x_2;
        handles.trajectory(2,i)=y_2;
        i=i+1;
        if i>400
            i=1;
        end
    elseif i~=1
        handles.trajectory=[];
        i=1;
    end
    [handles.t1,handles.t2,handles.w1,handles.w2]=db_pendulum(handles.t1,handles.t2,handles.w1,handles.w2,handles.mass1,handles.mass2,handles.len1,handles.len1,handles.G,.01);
    pend_plot(handles);
    set(handles.omega1,'string',num2str(handles.w1));
    set(handles.theta1,'string',num2str(handles.t1));
    set(handles.omega2,'string',num2str(handles.w2));
    set(handles.theta2,'string',num2str(handles.t2));
    guidata(hObject, handles);
end
set(hObject,'string','Start');
set(handles.omega1,'enable','on');
set(handles.theta1,'enable','on');
set(handles.omega2,'enable','on');
set(handles.theta2,'enable','on');
set(handles.length_rat,'enable','on');
set(handles.mass_rat,'enable','on');
set(handles.reset,'enable','on');
set(handles.gravity,'enable','on');
% --- Executes on slider movement.
function length_rat_Callback(hObject, eventdata, handles)
% hObject    handle to length_rat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.len1=get(hObject,'value');
handles.len2=1-handles.len1;
pend_plot(handles);
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function length_rat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to length_rat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function mass_rat_Callback(hObject, eventdata, handles)
% hObject    handle to mass_rat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mass1=get(hObject,'value');
handles.mass2=1-handles.mass1;
pend_plot(handles);
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function mass_rat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mass_rat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function omega1_Callback(hObject, eventdata, handles)
% hObject    handle to omega1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.w1=str2num(get(hObject,'String'));
energy(handles)
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of omega1 as text
%        str2double(get(hObject,'String')) returns contents of omega1 as a double


% --- Executes during object creation, after setting all properties.
function omega1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to omega1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theta1_Callback(hObject, eventdata, handles)
% hObject    handle to theta1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.t1=str2num(get(hObject,'String'));
pend_plot(handles);
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of theta1 as text
%        str2double(get(hObject,'String')) returns contents of theta1 as a double


% --- Executes during object creation, after setting all properties.
function theta1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function omega2_Callback(hObject, eventdata, handles)
% hObject    handle to omega2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.w1=str2num(get(hObject,'String'));
energy(handles)
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of omega2 as text
%        str2double(get(hObject,'String')) returns contents of omega2 as a
%        double


% --- Executes during object creation, after setting all properties.
function omega2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to omega2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theta2_Callback(hObject, eventdata, handles)
% hObject    handle to theta2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.t2=str2num(get(hObject,'String'));
pend_plot(handles);
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of theta2 as text
%        str2double(get(hObject,'String')) returns contents of theta2 as a double


% --- Executes during object creation, after setting all properties.
function theta2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in traj.
function traj_Callback(hObject, eventdata, handles)
% hObject    handle to traj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.trajectory=[];
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of traj


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.len1=.5;
handles.len2=.5;
handles.mass1=.5;
handles.mass2=.5;
handles.t1=0;
handles.t2=0;
handles.w1=0;
handles.w2=0;
handles.G=9.8;
set(handles.gr,'string','9.8')
set(handles.gravity,'value',9.8)
set(handles.omega1,'string','0')
set(handles.omega2,'string','0')
set(handles.theta1,'string','0')
set(handles.theta2,'string','0')
set(handles.mass_rat,'value',0.5)
set(handles.length_rat,'value',0.5)
pend_plot(handles);
guidata(hObject, handles);




function pend_plot(data);
[xx yy zz]=cylinder(0.1,20);
    x_1=data.len1*sin(data.t1);
    y_1=-data.len1*cos(data.t1) ;
    x_2=data.len2*sin(data.t2) + x_1;
    y_2=-data.len2*cos(data.t2) + y_1;
hold off
plot([0 x_1 x_2],[0 y_1 y_2],'k','linewidth',2);
hold on
if get(data.traj,'value') ==1
    plot(data.trajectory(1,:),data.trajectory(2,:),'r.','markersize',3)
end
fill(x_1+xx(1,:).*(data.mass1)^.3,y_1+yy(1,:).*(data.mass1)^.3,'b')
fill(x_2+xx(1,:).*(data.mass2)^.3,y_2+yy(1,:).*(data.mass2)^.3,'b')
set(gca,'xtick',[],'ytick',[])
axis(1.2*[-1 1 -1 1]);
box on;
drawnow
energy(data)



function [theta_1 theta_2 omega_1 omega_2]=db_pendulum(theta1, theta2, w1, w2, mass1, mass2, d1, d2, G, dt)
l1 = dt * (l(w1));
m1 = dt * (m(w2));
f1 = dt * (f(theta1, theta2, w1, w2, mass1, mass2, d1, d2, G));
g1 = dt * (g(theta1, theta2, w1, w2, mass1, mass2, d1, d2, G));

l2 = dt * (l(w1+(f1/2.)));
m2 = dt * (m(w2+(g1/2.)));
f2 = dt * (f(theta1+(l1/2.),theta2+(m1/2.),w1+(f1/2.),w2+(g1/2.), mass1, mass2, d1, d2, G));
g2 = dt * (g(theta1+(l1/2.),theta2+(m1/2.),w1+(f1/2.),w2+(g1/2.), mass1, mass2, d1, d2, G));

l3 = dt * (l(w1+(f2/2.)));
m3 = dt * (m(w2+(g2/2.)));
f3 = dt * (f(theta1+(l2/2.),theta2+(m2/2.),w1+(f2/2.),w2+(g2/2.), mass1, mass2, d1, d2, G));
g3 = dt * (g(theta1+(l2/2.),theta2+(m2/2.),w1+(f2/2.),w2+(g2/2.), mass1, mass2, d1, d2, G));

l4 = dt * (l(w1+f3));
m4 = dt * (m(w2+g3));
f4 = dt * (f(theta1+(l3),theta2+(m3),w1+(f3),w2+(g3), mass1, mass2, d1, d2, G));
g4 = dt * (g(theta1+(l3),theta2+(m3),w1+(f3),w2+(g3), mass1, mass2, d1, d2, G));
theta_1 = theta1 + (l1+(2.*l2)+(2.*l3)+l4)/6.;
theta_2 = theta2 + (m1+(2.*m2)+(2.*m3)+m4)/6.;
omega_1 = w1 + (f1+(2.*f2)+(2.*f3)+f4)/6.;
omega_2 = w2 + (g1+(2.*g2)+(2.*g3)+g4)/6.;

function f1=f(theta1, theta2, w1, w2, mass1, mass2, d1, d2, G)
f1= (  (( -d2*mass2*w2*w2*sin(theta1- theta2)/(d1*(mass1+mass2)) )  - ...
    (  mass2*w1*w1*cos(theta1-theta2)*sin(theta1- theta2)/(mass1+mass2) )  +...
    (  mass2*G*cos(theta1- theta2)*sin(theta2)/(d1*(mass1+mass2))  )  -...
    (  G*sin(theta1)/d1)...
    ) /  ( 1. - (mass2*cos(theta1- theta2)*cos(theta1- theta2)/(mass1+mass2)))  );
function g1=g(theta1, theta2 ,w1, w2, mass1, mass2, d1, d2, G)
g1=(( (-G*sin(theta2)/d2) + (mass2*w2*w2*cos(theta1-theta2)*sin(theta1-theta2)/(mass1+mass2)) +...
    (G*cos(theta1-theta2)*sin(theta1)/d2) + (d1*w1*w1*sin(theta1-theta2)/d2)  )  /  (1. - (mass2*cos(theta1-theta2)*cos(theta1-theta2)/(mass1+mass2))) );
function L1=l(w1)
L1=w1;
function M1=m(w2)
M1= w2;


% --- Executes on slider movement.
function gravity_Callback(hObject, eventdata, handles)
% hObject    handle to gravity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.G=get(hObject,'value');
set(handles.gr,'string',num2str(handles.G));
guidata(hObject, handles);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function gravity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gravity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function energy(data)
  Epot1=-data.mass1*data.G*data.len1*(cos(data.t1)-1.000) ;
  Epot2=-data.mass2*data.G*((cos(data.t1)-1.000)*data.len1+(cos(data.t2)-1.000)*data.len2);
  Ekin1=0.5*data.mass1*(data.w1*data.len1).^2;
  Ekin2=0.5*data.mass2*(   (data.w2*data.len2*cos(data.t2) +data.w1*data.len1*cos(data.t1)).^2 + (data.w2*data.len2*sin(data.t2) +data.w1*data.len1*sin(data.t1)).^2 );
  Etot=Ekin1+Ekin2+Epot1+Epot2;
  
  set(data.U1,'string',num2str(Epot1));
  set(data.U2,'string',num2str(Epot2));
  set(data.T1,'string',num2str(Ekin1));
  set(data.T2,'string',num2str(Ekin2));
  set(data.totE,'string',num2str(Etot));
