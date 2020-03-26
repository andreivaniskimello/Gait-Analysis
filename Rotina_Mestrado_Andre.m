% Routine to Drag Model.
% Calculate Drag Force during underwater walking from: anthropometric data and Lower Limb points position in time.
% By André Ivaniski Mello (andreivaniskimello@gmail.com)
% Last edit: 24.03.20

%% DATA ANALYSIS TO INSERT

% Winter residual analysis
% Phase portrait
% Fractal analysis?
% Joints Angles


%%
%%Orientations%%

%1º) Insert SubjectName_Depth__Speed_StrideNumber
%2º) Insert files path of 'Anthropometric' and 'Linear' datas
%3º) Insert Output file path
%4º) Insert TD1, TO, TD2 frames
%5º) Run routine


% THIS ROUTINE HAS 12 SECTIONS

% SECTION 01: INFORMATION FOR ROUTINE AND INPUT DATA IMPORTATION
% SECTION 02: DATA INFORMATION PREPARATION 
% SECTION 03: SPATIOTEMPORAL STRIDE PARAMETERS
% SECTION 04: ANGULAR STRIDE PARAMETERS
% SECTION 05: DRAG FORCE MODEL DATA PREPARATION
% SECTION 06: XIPHOID TRUNK MODEL
% SECTION 07: UMBILICAL TRUNK MODEL
% SECTION 08: UPPER LEG (UL) MODEL
% SECTION 09: LOWER LEG (LL) MODEL
% SECTION 10: FOOT (FT) MODEL
% SECTION 11: GROUND REACTION FORCES ESTIMATION
% SECTION 12: ALGORITHM FOR AUTOMATIC TOUCH DOWN AND TAKE OFF EVENTS %
% SECTION 13: EXPORT DATA
% SECTION 14: GRAPHICS AND GIFS

clear all
close all
clc

%% SECTION 01
% INFORMATION FOR ROUTINE AND INPUT DATA IMPORTATION %

File_Name =['Gabriel_Xifoide_0.4_P3'];        %insert 'SubjectName_Depth_Speed_StrideNumber' here, in order to export the final data with the Subject Name_Depth 

%Open Anthropometric Values
Anthropometric_file = fopen('C:\Users\andre\Documents\Mestrado\Projetos Pesquisa\Projeto Marcha Água\Coletas\Data\Cinemetria\Antropometria\Gabriel_Antropometria.txt');
data = textscan(Anthropometric_file,'%f%s');
fclose(Anthropometric_file);
Anthropometric = data{1,1};

%Open Walk Linear Data
Walk_Linear_file = ('C:\Users\andre\Documents\Mestrado\Projetos Pesquisa\Projeto Marcha Água\Coletas\Data\Cinemetria\Data Out\Gabriel\Xifoide\Gabriel_Xifoide_0.4_P3_Positions.txt');
delimiterIn = '\t';
Linear_Data = importdata(Walk_Linear_file,delimiterIn);
Linear_Data = Linear_Data.';

%Insert Export (Output) File Path here
Export_File_Path = ['C:\Users\andre\Documents\Mestrado\Projetos Pesquisa\Projeto Marcha Água\Coletas\Data\Cinemetria\Data Out\Gabriel\Xifoide'];     

% Insert TD and TO frames 
    % TD: First frame that the foot touch the ground
    % TO: First frame that the foot is out of the ground

TD1 = 25;
TO = 91;                              
TD2 = 138;

TO = TO - TD1;     % This is adjusting TO frame number considering TD1 as Frame 1

%% SECTION 02
% DATA INFORMATION PREPARATION %

    % ANTHROPOMETRIC DATA
    
Anthropometric = Anthropometric/100;       %Convert anthropometric values in centimeters (cm) to meters (m)

%Anthropometric parameters labeling

Stature = Anthropometric(1);
ShoulderH = Anthropometric(2);
SubsternalH = Anthropometric(3);
UmbilicalH = Anthropometric(4);
SittingH = Anthropometric(5);
TrochantericH = Anthropometric(6);
TibialH = Anthropometric(7);
SphrH = Anthropometric(8);
FootL = Anthropometric(9);
XiphoidC = Anthropometric(10);
UmbilicalC = Anthropometric(11);
HipB = Anthropometric(12);
ThighC = Anthropometric(13);
KneeC = Anthropometric(14);
AnkleC = Anthropometric(15);


    %POSITION DATA

%%Linear data file in the format (n,10)
% n = frames
% [Toe Heel Ankle Knee Trochanter Trunk]
% 02 columns per Point: x y
% x: horizontal; y: vertical


Linear_Data = Linear_Data(TD1:TD2-1,:);

Frames = size(Linear_Data,1);                    % Calculate the total frames number
Toe_Raw = Linear_Data(:,1:2);                     
Heel_Raw = Linear_Data(:,3:4);
Ankle_Raw = Linear_Data(:,5:6);
Knee_Raw = Linear_Data(:,7:8);
Trochanter_Raw = Linear_Data(:,9:10);
Trunk_Point_Raw = Linear_Data(:,11:12);

Frames_Vector = [1:1:Frames];
Frames_Vector = Frames_Vector.';

%%% LOWPASS BUTTERWORTH FILTER %%%

fsample = 60;
dt=1/fsample;
fcut=6;
order=2;

Toe = matfiltfilt(dt, fcut, order, Toe_Raw);
Heel = matfiltfilt(dt, fcut, order, Heel_Raw);
Ankle = matfiltfilt(dt, fcut, order, Ankle_Raw);
Knee = matfiltfilt(dt, fcut, order, Knee_Raw);
Trochanter = matfiltfilt(dt, fcut, order, Trochanter_Raw);
Trunk_Point = matfiltfilt(dt, fcut, order, Trunk_Point_Raw);

% Points Axis Labeling
Toe_x = Toe(:,1);
Toe_y = Toe(:,2);
Heel_x = Heel(:,1);
Heel_y = Heel(:,2);
Ankle_x = Ankle(:,1);
Ankle_y = Ankle(:,2);
Knee_x = Knee(:,1);
Knee_y = Knee(:,2);
Trochanter_x = Trochanter(:,1);
Trochanter_y = Trochanter(:,2);
Trunk_Point_x = Trunk_Point(:,1);
Trunk_Point_y = Trunk_Point(:,2);

%% SECTION 03
% SPATIOTEMPORAL STRIDE PARAMETERS %

Time_Total = Frames/fsample;                                     %Time = total frames / fsample
Time_Vector = [0:dt:((TD2-TD1)/fsample)];                        %Create a Time Column (Vector)
Time_Stride_Percentage = Time_Vector/(max(Time_Vector))*100;     %Create a Time Vector (row) (1xn) in Percentage of the Stride
Time_Stride_Percentage = Time_Stride_Percentage.';               %Rotate Time Vector (column) in Percentage of the Stride to (nx1)
[nlTimeStride, ~]=size(Time_Stride_Percentage);
TO_in_Percentage_Decimal = TO/(nlTimeStride)*100;
TO_in_Percentage = fix(TO_in_Percentage_Decimal);                 %Calculate TO in Percentage of the Stride

Stride_Lenght = (Heel_x(end) - Heel_x(1));                        %Final position - Initial Position
Speed_Mean = (Stride_Lenght/Time_Total);                          %Speed_Mean: is the mean speed value over the entire walk.
Stride_Frequency = Speed_Mean/Stride_Lenght;

Stride_Duration = (TD2 - TD1)*dt;
Contact_Phase_Duration = (TO - TD1)*dt;
Duty_Factor = (Contact_Phase_Duration/Stride_Duration)*100;

%% SECTION 04
% ANGULAR STRIDE PARAMETERS %

% SEGMENTS VECTORS DETERMINATION
Foot = Heel - Toe;
Shank = Knee - Ankle;
Thigh = Trochanter - Knee;
Trunk = Trunk_Point - Trochanter;

%The loops below corrects the SEGMENTS VECTOR to keep them only at positive
%cartesian quadrant. This is to prevent any disruption at the ANGLES
%curves.

% Foot
for i=1:Frames    
    if (min(Foot(:,1))) <= 0
        Foot_Corrected(i,1) = (Foot(i,1)-(min(Foot(:,1))))+ (abs(min(Foot(:,1))));
        Foot_Corrected(i,2) =  Foot_Corrected(i,1) * abs((Foot(i,2)./Foot(i,1)));
        
    elseif (min(Foot(:,2))) <= 0
        Foot_Corrected(i,2) = (Foot(i,2)-(min(Foot(:,2))))+ (abs(min(Foot(:,2))));
        Foot_Corrected(i,1) = Foot_Corrected(i,2) * abs((Foot(i,2)./Foot(i,1)));
        
    else
        Foot_Corrected(i,:) = Foot(i,:);
    end
end


% Shank
for i=1:Frames    
    if (min(Shank(:,1))) <= 0
        Shank_Corrected(i,1) = (Shank(i,1)-(min(Shank(:,1))))+ (abs(min(Shank(:,1))));
        Shank_Corrected(i,2) = Shank_Corrected(i,1) * abs((Shank(i,2)./Shank(i,1)));
    
    elseif (min(Shank(:,2))) <= 0
        Shank_Corrected(i,2) = (Shank(i,2)-(min(Shank(:,2))))+ (abs(min(Shank(:,2))));
        Shank_Corrected(i,1) = Shank_Corrected(i,2) * abs(Shank(i,2)./Shank(i,1));
    
    else
        Shank_Corrected(i,:) = Shank(i,:);
    end
end


% Thigh
for i=1:Frames    
    if (min(Thigh(:,1))) <= 0
        Thigh_Corrected(i,1) = (Thigh(i,1)-(min(Thigh(:,1))))+ (abs(min(Thigh(:,1))));
        Thigh_Corrected(i,2) =  Thigh_Corrected(i,1) * abs((Thigh(i,2)./Thigh(i,1)));
        
    elseif (min(Thigh(:,2))) <= 0
        Thigh_Corrected(i,2) = (Thigh(i,2)-(min(Thigh(:,2))))+ (abs(min(Thigh(:,2))));
        Thigh_Corrected(i,1) = Thigh_Corrected(i,2) * abs((Thigh(i,2)./Thigh(i,1)));
        
    else
        Thigh_Corrected(i,:) = Thigh(i,:);    
    end
end


% Trunk
for i=1:Frames    
    if (min(Trunk(:,1))) <= 0
        Trunk_Corrected(i,1) = (Trunk(i,1)-(min(Trunk(:,1))))+ (abs(min(Trunk(:,1))));
        Trunk_Corrected(i,2) =  Trunk_Corrected(i,1) * abs((Trunk(i,2)./Trunk(i,1)));
    
    elseif (min(Trunk(:,2))) <= 0
        Trunk_Corrected(i,2) = (Trunk(i,2)-(min(Trunk(:,2))))+ (abs(min(Trunk(:,2))));
        Trunk_Corrected(i,1) = Trunk_Corrected(i,2) * abs((Trunk(i,2)./Trunk(i,1)));

    else
        Trunk_Corrected(i,:) = Trunk(i,:);    
    end
end


    % SEGMENT'S ANGULAR POSITION

Foot_Angle = atand(Foot_Corrected(:,2)./Foot_Corrected(:,1));
Shank_Angle = atand(Shank_Corrected(:,2)./Shank_Corrected(:,1));
Thigh_Angle = atand(Thigh_Corrected(:,2)./Thigh_Corrected(:,1));
Trunk_Angle = atand(Trunk_Corrected(:,2)./Trunk_Corrected(:,1));

    % SEGMENT'S ANGULAR RANGE OF MOTION (ROM)
    
Foot_Angle_ROM = max(Foot_Angle)- min(Foot_Angle);
Shank_Angle_ROM = max(Shank_Angle)- min(Shank_Angle);
Thigh_Angle_ROM = max(Thigh_Angle)- min(Thigh_Angle);
Trunk_Angle_ROM = max(Trunk_Angle)- min(Trunk_Angle);


    % SEGMENT'S ANGULAR SPEED

% Foot
for i = 2:Frames-1
    Foot_Angular_Speed (i) = (Foot_Angle(i+1) - Foot_Angle(i-1))/(2*dt);
end
Foot_Angular_Speed = interpft(Foot_Angular_Speed,Frames);
Foot_Angular_Speed = Foot_Angular_Speed.';
Foot_Angular_Speed_Max = max(Foot_Angular_Speed);

% Shank
for i = 2:Frames-1
    Shank_Angular_Speed (i) = (Shank_Angle(i+1) - Shank_Angle(i-1))/(2*dt);
end
Shank_Angular_Speed = interpft(Shank_Angular_Speed,Frames);
Shank_Angular_Speed = Shank_Angular_Speed.';
Shank_Angular_Speed_Max = max(Shank_Angular_Speed);

% Thigh
for i = 2:Frames-1
    Thigh_Angular_Speed (i) = (Thigh_Angle(i+1) - Thigh_Angle(i-1))/(2*dt);
end
Thigh_Angular_Speed = interpft(Thigh_Angular_Speed,Frames);
Thigh_Angular_Speed = Thigh_Angular_Speed.';
Thigh_Angular_Speed_Max = max(Thigh_Angular_Speed);

% Trunk
for i = 2:Frames-1
    Trunk_Angular_Speed (i) = (Trunk_Angle(i+1) - Trunk_Angle(i-1))/(2*dt);
end
Trunk_Angular_Speed = interpft(Trunk_Angular_Speed,Frames);
Trunk_Angular_Speed = Trunk_Angular_Speed.';
Trunk_Angular_Speed_Max = max(Trunk_Angular_Speed);



    % JOINT'S ANGULAR POSITION

Ankle_Joint_Angle = Shank_Angle - Foot_Angle;
Knee_Joint_Angle = Thigh_Angle - Shank_Angle;
Hip_Joint_Angle = Trunk_Angle - Thigh_Angle;


    % JOINT'S ANGULAR RANGE OF MOTION (ROM)
    
Ankle_Joint_Angle_ROM = max(Ankle_Joint_Angle)- min(Ankle_Joint_Angle);
Knee_Joint_Angle_ROM = max(Knee_Joint_Angle)- min(Knee_Joint_Angle);
Hip_Joint_Angle_ROM = max(Hip_Joint_Angle)- min(Hip_Joint_Angle);
   

    % JOINT'S ANGULAR SPEED

% Ankle Joint
for i = 2:Frames-1
    Ankle_Joint_Angular_Speed (i) = (Ankle_Joint_Angle(i+1) - Ankle_Joint_Angle(i-1))/(2*dt);
end
Ankle_Joint_Angular_Speed = interpft(Ankle_Joint_Angular_Speed, Frames);
Ankle_Joint_Angular_Speed = Ankle_Joint_Angular_Speed.';
Ankle_Joint_Angular_Speed_Max = max(Ankle_Joint_Angular_Speed);

% Knee Joint
for i = 2:Frames-1
    Knee_Joint_Angular_Speed (i) = (Knee_Joint_Angle(i+1) - Knee_Joint_Angle(i-1))/(2*dt);
end
Knee_Joint_Angular_Speed = interpft(Knee_Joint_Angular_Speed, Frames);
Knee_Joint_Angular_Speed = Knee_Joint_Angular_Speed.';
Knee_Joint_Angular_Speed_Max = max(Knee_Joint_Angular_Speed);

% Hip Joint
for i = 2:Frames-1
    Hip_Joint_Angular_Speed (i) = (Hip_Joint_Angle(i+1) - Hip_Joint_Angle(i-1))/(2*dt);
end
Hip_Joint_Angular_Speed = interpft(Hip_Joint_Angular_Speed, Frames);
Hip_Joint_Angular_Speed = Hip_Joint_Angular_Speed.';
Hip_Joint_Angular_Speed_Max = max(Hip_Joint_Angular_Speed);



    % PHASE PORTRAIT

% Foot
Foot_Angle_Normalized = Foot_Angle/max(Foot_Angle(:,1));
Foot_Angular_Speed_Normalized  = Foot_Angular_Speed/max(Foot_Angular_Speed(:,1));
Foot_Phase_Angle = atand(Foot_Angular_Speed_Normalized(:,1)./Foot_Angle_Normalized(:,1));      

% Shank
Shank_Angle_Normalized = Shank_Angle/max(Shank_Angle(:,1));
Shank_Angular_Speed_Normalized  = Shank_Angular_Speed/max(Shank_Angular_Speed(:,1));
Shank_Phase_Angle = atand(Shank_Angular_Speed_Normalized(:,1)./Shank_Angle_Normalized(:,1)); 

% Thigh
Thigh_Angle_Normalized = Thigh_Angle/max(Thigh_Angle(:,1));
Thigh_Angular_Speed_Normalized  = Thigh_Angular_Speed/max(Thigh_Angular_Speed(:,1));
Thigh_Phase_Angle = atand(Thigh_Angular_Speed_Normalized(:,1)./Thigh_Angle_Normalized(:,1)); 


    % RELATIVE PHASE
Relative_Phase_Shank_Thigh = Shank_Phase_Angle - Thigh_Phase_Angle ;


%% SECTION 05
% DRAG FORCE MODEL DATA PREPARATION %

%Calculate instantenous Velocity Vector in x Axis of each Anatomical point by Finite Difference Technique Method (Winter, 2009)%

%Trunk Velocity
for i = 2:(length(Trunk_Point)-1)
    Velocity_Trunk(i) = (Trunk_Point_x(i+1) - Trunk_Point_x(i-1))/(2*dt);
end

%Throcanter Velocity
for i = 2:(length(Trochanter)-1)
    Velocity_Throcanter(i) = (Trochanter_x(i+1)- Trochanter_x(i-1))/(2*dt);
end

%Knee Velocity
for i = 2:(length(Knee)-1)
    Velocity_Knee(i) = (Knee_x(i+1)- Knee_x(i-1))/(2*dt);
end

%Ankle Velocity
for i = 2:(length(Ankle_x)-1)
    Velocity_Ankle(i) = (Ankle_x(i+1)- Ankle_x(i-1))/(2*dt);
end

%Heel Velocity
for i = 2:(length(Heel_x)-1)
    Velocity_Heel(i) = (Heel_x(i+1)- Heel_x(i))/(2*dt);
end
        
%Toe Velocity
for i = 2:(length(Toe_x)-1)
    Velocity_Toe(i) = (Toe_x(i+1)- Toe_x(i-1))/(2*dt);
end

 
V_Trunk_Mean = mean(Velocity_Trunk);             %Just asking for Mean Velocity for overall analysis.
V_Throcanter_Mean = mean(Velocity_Throcanter);
V_Knee_Mean = mean(Velocity_Knee);
V_Ankle_Mean = mean(Velocity_Ankle);
V_Toe_Mean = mean(Velocity_Toe);
V_Heel_Mean = mean(Velocity_Heel);


%%% GENERAL FUNCTIONS FOR DRAG FORCE %%%

%%dA(z) = 2*((Ad-Ap)z + Ap*L)/L))*dz
%%v(z)=((z*Vd)+(L-z)*Vp)/L
%%dAP(z) = d(z).sen(angle)   , dAP: Projected Area at Vertical Axis
%%(perpendicular to Velocity Vector)

%% SECTION 06
% XIPHOID TRUNK MODEL %

%Ap = HipBreadth/2
%Ad = XiphoidC/2(pi)
%L = SittingH - (Stature - SubsternalH)

Ap_Xiph = HipB/2;
Ad_Xiph = XiphoidC/(2*pi);
L_Xiph = SittingH - (Stature - SubsternalH);

%Xiphoid Trunk length differential
dL_Xiph = 0:0.001:max(L_Xiph);
 
%%dA(z) = 2*((Ad-Ap)z + Ap*L)/L))*dz
% z = (0+dUL_L(i))

%dA_Xiph = dA(z)

for i=1:length(dL_Xiph)
    dA_Xiph(i) =(((Ad_Xiph-Ap_Xiph)*(0+dL_Xiph(i))+(Ap_Xiph*L_Xiph))/L_Xiph)*0.001*2;
end

%dAP_Xiph [frames dA_Xiph parts] = dAP(z)
for i=1:length(Trunk_Angle);
    
    for j=1:length(dA_Xiph);
        if Trunk_Angle(i)<=90;
           dAP_Xiph(i,j)=dA_Xiph(j)*sind(Trunk_Angle(i));
        else
            dAP_Xiph(i,j)=dA_Xiph(j)*sind(180-Trunk_Angle(i));
        end
    end
end

%%v(z)=((z*Vd)+(L-z)*Vp)/L

%Xiph Velocity
%Vp: Velocity_Throcanter
%Vd: Velocity_Xiphoid
%V_Xiph = v(z)

V_Xiph = zeros(length(Velocity_Throcanter), length(dL_Xiph));   %preallocation: to optimize routine running time


%Nested For Loop. First For Loop to time vector. Second For Loop to Segment
%dL(z).

%V_Xiph [frames dL_Xiph]

for j=1:1:length(Velocity_Throcanter)
    
    for q=1:length(dL_Xiph)
        V_Xiph (j,q) = (((0+dL_Xiph(q))*Velocity_Trunk(j))+((L_Xiph-(0+dL_Xiph(q)))*Velocity_Throcanter(j)))/L_Xiph;
    end
end

%Cd: drag coefficient;
%W: water density (values from Orselli & Duarte, 2011)
Cd = 1;
W = 1000;

for j=1:size(V_Xiph,1)
    for q=1:size(V_Xiph,2)
        dF_Xiph(j,q) = Cd*W*0.5*(V_Xiph(j,q).^2)*(dAP_Xiph(j,q));
    end
end

F_Xiph = (-1)*sum(dF_Xiph,2);

F_Xiph_Filtered = matfiltfilt(dt, fcut, order, F_Xiph);


%Dividing Drag Force in Contact and Swing Phases

F_Xiph_Contact = F_Xiph_Filtered(1:TO-1);
F_Xiph_Swing = F_Xiph_Filtered(TO:length(F_Xiph_Filtered));

F_Xiph_Contact_Max = max(F_Xiph_Contact);
F_Xiph_Contact_Min = min(F_Xiph_Contact);
F_Xiph_Contact_Mean = mean(F_Xiph_Contact);

F_Xiph_Swing_Max = max(F_Xiph_Swing);
F_Xiph_Swing_Min = min(F_Xiph_Swing);
F_Xiph_Swing_Mean = mean(F_Xiph_Swing);

%Calulating RMS of Drag Force during Contact and Swing Phases

RMS_F_Xiph_Contact = rms(F_Xiph_Contact);
RMS_F_Xiph_Swing = rms(F_Xiph_Swing);

figure ('Name', 'Drag Force Trunk')
title ('Drag Force Trunk');
xlabel('Frames (n)');
ylabel ('Drag Force (N)');
hold on
plot(F_Xiph,'r');
hold on;
plot(F_Xiph_Filtered,'k','LineWidth',2.0);
legend ('Raw', 'Filtered LowPass Butter');
hold on

%% SECTION 07
% UMBILICAL TRUNK MODEL %

%Ap = HipBreadth/2
%Ad = UmbilicalC/2(pi)
%L = SittingH - (Stature - UmbilicalH)

Ap_Umbi = HipB/2;
Ad_Umbi = UmbilicalC/(2*pi);
L_Umbi = SittingH - (Stature - UmbilicalH);

%Umbilical Trunk length differential
dL_Umbi = 0:0.001:max(L_Umbi);

%%dA(z) = 2*((Ad-Ap)z + Ap*L)/L))*dz
% z = (0+dUL_L(i))

%dA_Umbi = dA(z)

for i=1:length(dL_Umbi)
    dA_Umbi(i) =(((Ad_Umbi-Ap_Umbi)*(0+dL_Umbi(i))+(Ap_Umbi*L_Umbi))/L_Umbi)*0.001*2;
end

%dAP_Umbi [frames dA_Umbi parts] = dAP(z)
for i=1:length(Trunk_Angle);
    for j=1:length(dA_Umbi);
        if Trunk_Angle(i)<=90;
            dAP_Umbi(i,j)=dA_Umbi(j)*sind(Trunk_Angle(i));
        else
            dAP_Umbi(i,j)=dA_Umbi(j)*sind(180-Trunk_Angle(i));
        end
    end
end

%%v(z)=((z*Vd)+(L-z)*Vp)/L

%Umbi Velocity
%Vp: Velocity_Throcanter
%Vd: Velocity_Umbilical
%V_Umbi = v(z)

V_Umbi = zeros(length(Velocity_Throcanter), length(dL_Umbi));   %preallocation: to optimize routine running time


%Nested For Loop. First For Loop to time vector. Second For Loop to Segment
%dL(z).

%V_Xiph [frames dL_Xiph]

for j=1:1:length(Velocity_Throcanter)
    for q=1:length(dL_Umbi)
        V_Umbi (j,q) = (((0+dL_Umbi(q))*Velocity_Trunk(j))+((L_Umbi-(0+dL_Umbi(q)))*Velocity_Throcanter(j)))/L_Umbi;
    end
end

%Cd: drag coefficient;
%W: water density (values from Orselli & Duarte, 2011)
Cd = 1;
W = 1000;

for j=1:size(V_Umbi,1)
    for q=1:size(V_Umbi,2)
        dF_Umbi(j,q) = Cd*W*0.5*(V_Umbi(j,q).^2)*(dAP_Umbi(j,q));
    end
end

F_Umbi = (-1)*sum(dF_Umbi,2);

F_Umbi_Filtered = matfiltfilt(dt, fcut, order, F_Umbi);


%Dividing Drag Force in Contact and Swing Phases

F_Umbi_Contact = F_Umbi_Filtered(1:TO-1);
F_Umbi_Swing = F_Umbi_Filtered(TO:length(F_Umbi_Filtered));

F_Umbi_Contact_Max = max(F_Umbi_Contact);
F_Umbi_Contact_Min = min(F_Umbi_Contact);
F_Umbi_Contact_Mean = mean(F_Umbi_Contact);

F_Umbi_Swing_Max = max(F_Umbi_Swing);
F_Umbi_Swing_Min = min(F_Umbi_Swing);
F_Umbi_Swing_Mean = mean(F_Umbi_Swing);


%Calulating RMS of Drag Force during Contact and Swing Phases

RMS_F_Umbi_Contact = rms(F_Umbi_Contact);
RMS_F_Umbi_Swing = rms(F_Umbi_Swing);

%% SECTION 08
% UPPER LEG (UL) MODEL %

%Ap = ThighC/2(pi)
%Ad = KneeC/2(pi)
%L = Stature - SittingH - TibialH
%DELSH = SittinhH - Stature - TrochantericH
%DELSH: distance from upper segment end to joint center.

%Upper Leg Size Modeling (Cylinder)
Ap_UL = ThighC/(2*pi);
Ad_UL = KneeC/(2*pi);
L_UL = Stature - SittingH - TibialH;

%UL length differential
dL_UL=0:0.001:max(L_UL);


%%dA(z) = 2*((Ad-Ap)z + Ap*L)/L))*dz
% z = (0+dUL_L(i))

%dA_UL = dA(z)

for i=1:length(dL_UL)
    dA_UL(i) =(((Ap_UL-Ad_UL)*(0+dL_UL(i))+(Ad_UL*L_UL))/L_UL)*0.001*2;
end

%dAP_UL [frames dA_ULparts] = dAP(z)
for i=1:length(Thigh_Angle);
    
    for j=1:length(dA_UL);
        if Thigh_Angle(i)<=90;
            dAP_UL(i,j)=dA_UL(j)*sind(Thigh_Angle(i));
        else
            dAP_UL(i,j)=dA_UL(j)*sind(180-Thigh_Angle(i));
        end
    end
end


%%v(z)=((z*Vd)+(L-z)*Vp)/L

%UL Velocity
%Vp: Velocity_Iliaca
%Vd: Velocity_Joelho
%V_UL = v(z)

V_UL = zeros(length(Velocity_Knee), length(dL_UL));   %preallocation: to optimize routine running time


%Nested For Loop. First For Loop to time vector. Second For Loop to Segment
%dL(z).

%V_UL [frames dL_UL]

for j=1:1:length(Velocity_Knee)
    for q=1:length(dL_UL)
        V_UL (j,q) = (((0+dL_UL(q))*Velocity_Throcanter(j))+((L_UL-(0+dL_UL(q)))*Velocity_Knee(j)))/L_UL;
    end
end


%Cd: drag coefficient;
%W: water density (values from Orselli & Duarte, 2011)
Cd = 1;
W = 1000;

for j=1:size(V_UL,1)
    for q=1:size(V_UL,2)
        dF_UL(j,q) = Cd*W*0.5*(V_UL(j,q).^2)*(dAP_UL(j,q));
    end
end

F_UL = (-1)*sum(dF_UL,2);

F_UL_Filtered = matfiltfilt(dt, fcut, order, F_UL);


%Dividing Drag Force in Contact and Swing Phases

F_UL_Contact = F_UL_Filtered(1:TO-1);
F_UL_Swing = F_UL_Filtered(TO:length(F_UL_Filtered));

F_UL_Contact_Max = max(F_UL_Contact);
F_UL_Contact_Min = min(F_UL_Contact);
F_UL_Contact_Mean = mean(F_UL_Contact);

F_UL_Swing_Max = max(F_UL_Swing);
F_UL_Swing_Min = min(F_UL_Swing);
F_UL_Swing_Mean = mean(F_UL_Swing);


%Calulating RMS of Drag Force during Contact and Balance Phases

RMS_F_UL_Contact = rms(F_UL_Contact);
RMS_F_UL_Swing = rms(F_UL_Swing);


figure ('Name', 'Drag Force Upper Leg')
title ('Drag Force Upper Leg');
xlabel('Frames (n)');
ylabel ('Drag Force (N)');
hold on
plot(F_UL,'r');
hold on;
plot(F_UL_Filtered,'k','LineWidth',2.0);
legend ('Raw', 'Filtered LowPass Butter');
hold on

%% SECTION 09
% LOWER LEG (LL) MODEL %

%Ap = KneeC/2(pi)
%Ad = AnkleC/2(pi)
%SL = TibialH - SphrH

%Lower Leg Size Modeling (Cylinder)
Ap_LL = KneeC/(2*pi);
Ad_LL = AnkleC/(2*pi);
L_LL = TibialH - SphrH;

%LL length differential
dL_LL=0:0.001:max(L_LL);


%%dA(z) = 2*((Ad-Ap)z + Ap*L)/L))*dz
% z = (0+dL_LL(i))

%dA_LL = dA(z)

for i=1:length(dL_LL)
    dA_LL(i) =(((Ap_LL-Ad_LL)*(0+dL_LL(i))+(Ad_LL*L_LL))/L_LL)*0.001*2;
end

%dAP_LL [frames dA_LLparts] = dAP(z)
for i=1:length(Shank_Angle);
    
    for j=1:length(dA_LL);
        if Shank_Angle(i)<=90;
            dAP_LL(i,j)=dA_LL(j)*sind(Shank_Angle(i));
        else
            dAP_LL(i,j)=dA_LL(j)*sind(180-Shank_Angle(i));
        end
    end
end


%%v(z)=((z*Vd)+(L-z)*Vp)/L

%LL Velocity
%Vp: Velocity_Joelho
%Vd: Velocity_Maleolo
%V_UL = v(z)

V_LL = zeros(length(Velocity_Knee), length(dL_LL));   %preallocation: to optimize routine running time


%Nested For Loop. First For Loop to time vector. Second For Loop to Segment
%dL(z).

%V_UL [frames dL_LL]

for j=1:1:length(Velocity_Knee)
    
    for q=1:length(dL_LL)
        V_LL (j,q) = (((0+dL_LL(q))*Velocity_Ankle(j))+((L_LL-(0+dL_LL(q)))*Velocity_Knee(j)))/L_LL;
    end
end

%Cd: drag coefficient;
%W: water density (values from Orselli & Duarte, 2011)
Cd = 1;
W = 1000;

for j=1:size(V_LL,1)
    for q=1:size(V_LL,2)
        dF_LL(j,q) = Cd*W*0.5*(V_LL(j,q).^2)*(dAP_LL(j,q));
    end
end

F_LL = (-1)*sum(dF_LL,2);


%y = lowpass(x,fpass,fs) specifies that x is sampled at a rate of fs hertz. fpass is the passband frequency of the filter in hertz.

dt=1/60;
fcut=8;
order=2;

F_LL_Filtered = matfiltfilt(dt, fcut, order, F_LL);


%Dividing Drag Force in Contact and Swing Phases

F_LL_Contact = F_LL_Filtered(1:TO-1);
F_LL_Swing = F_LL_Filtered(TO:length(F_LL_Filtered));

F_LL_Contact_Max = max(F_LL_Contact);
F_LL_Contact_Min = min(F_LL_Contact);
F_LL_Contact_Mean = mean(F_LL_Contact);

F_LL_Swing_Max = max(F_LL_Swing);
F_LL_Swing_Min = min(F_LL_Swing);
F_LL_Swing_Mean = mean(F_LL_Swing);


%Calulating RMS of Drag Force during Contact and Swing Phases

RMS_F_LL_Contact = rms(F_LL_Contact);
RMS_F_LL_Swing = rms(F_LL_Swing);

figure ('Name', 'Drag Force Lower Leg')
title ('Drag Force Lower Leg');
xlabel('Frames (n)');
ylabel ('Drag Force (N)');
hold on
plot(F_LL,'r');
hold on;
plot(F_LL_Filtered,'k','LineWidth',2.0);
legend ('Raw', 'Filtered LowPass Butter');
hold off;

%% SECTION 10
% FOOT (FT) MODEL %

%Ap = SphrH*0.5
%Ad = ((1.287*L)-L)/((1.287*L)/Ap_Ft)           %Calculate from trigonometry (Triangles Similarity), CG from Foot is at Ap/3 from axis, and at 0.429L from Proximal Joint.
%SL = FootL

%Foot Size Modeling (Cylinder)
Ap_FT = SphrH*0.5;
L_FT = FootL;
Ad_FT = ((1.287*L_FT)-L_FT)/((1.287*L_FT)/Ap_FT);


%FT length differential
dL_FT=0:0.001:max(L_FT);


%%dA(z) = 2*((Ad-Ap)z + Ap*L)/L))*dz
% z = (0+dL_FT(i))

%dA_FT = dA(z)

for i=1:length(dL_FT)
    dA_FT(i) =(((Ap_FT-Ad_FT)*(0+dL_FT(i))+(Ad_FT*L_FT))/L_FT)*0.001*2;
end

%dAP_FT [frames dA_FTparts] = dAP(z)
for i=1:length(Foot_Angle);
    
    for j=1:length(dA_FT);
        if Foot_Angle(i)<=90;
            dAP_FT(i,j)=dA_FT(j)*sind(Foot_Angle(i));
        else
            dAP_FT(i,j)=dA_FT(j)*sind(180-Foot_Angle(i));
        end
    end
end


%%v(z)=((z*Vd)+(L-z)*Vp)/L

%FT Velocity
%Vp: Velocity_Maleolo
%Vd: Velocity_Meta
%V_UL = v(z)

V_FT = zeros(length(Velocity_Ankle), length(dL_FT));   %preallocation: to optimize routine running time


%Nested For Loop. First For Loop to time vector. Second For Loop to Segment
%dL(z).

%V_FT [frames dL_FT]

for j=1:1:length(Velocity_Ankle)
    
    for q=1:length(dL_FT)
        V_FT (j,q) = (((0+dL_FT(q))*Velocity_Toe(j))+((L_FT-(0+dL_FT(q)))*Velocity_Ankle(j)))/L_FT;
    end
end

%Cd: drag coefficient;
%W: water density (values from Orselli & Duarte, 2011)
Cd = 1;
W = 1000;

for j=1:size(V_FT,1)
    for q=1:size(V_FT,2)
        dF_FT(j,q) = Cd*W*0.5*(V_FT(j,q).^2)*(dAP_FT(j,q));
    end
end

F_FT = (-1)*sum(dF_FT,2);

%y = lowpass(x,fpass,fs) specifies that x is sampled at a rate of fs hertz. fpass is the passband frequency of the filter in hertz.

dt=1/60;
fcut=8;
order=2;

F_FT_Filtered = matfiltfilt(dt, fcut, order, F_FT);

%Dividing Drag Force in Contact and Swing Phases

F_FT_Contact = F_FT_Filtered(1:TO-1);
F_FT_Swing = F_FT_Filtered(TO:length(F_FT_Filtered));

F_FT_Contact_Max = max(F_FT_Contact);
F_FT_Contact_Min = min(F_FT_Contact);
F_FT_Contact_Mean = mean(F_FT_Contact);


F_FT_Swing_Max = max(F_FT_Swing);
F_FT_Swing_Min = min(F_FT_Swing);
F_FT_Swing_Mean = mean(F_FT_Swing);

%Calulating RMS of Drag Force during Contact and Balance Phases

RMS_F_FT_Contact = rms(F_FT_Contact);
RMS_F_FT_Swing = rms(F_FT_Swing);

figure ('Name', 'Drag Force Foot')
title ('Drag Force Foot');
xlabel('Frames (n)');
ylabel ('Drag Force (N)');
hold on
plot(F_FT,'r');
hold on;
plot(F_FT_Filtered,'k','LineWidth',2.0);
legend ('Raw', 'Filtered LowPass Butter');
hold off;

%%

figure('Name','Segments Angles','ScreenSize')
title('Angles');
hold on
plot(Foot_Angle(:,1),'r','LineWidth',2.5)
hold on
plot (Shank_Angle(:,1),'b', 'LineWidth',2.5)
hold on
plot (Thigh_Angle(:,1), 'k','LineWidth',2.5)
hold on
plot (Trunk_Angle(:,1), 'g','LineWidth',2.5)
hold on
plot([TO TO],[0 100],'-m', 'LineWidth',2.5)
hold on
plot([0 120],[90 90],'-k', 'LineWidth',2.5)
legend ('Foot', 'Shank', 'Thigh','Trunk', 'TO','Vertical')

figure('Name','Segments Angular Speed')
title('Angular Speed');
hold on
plot(Foot_Angular_Speed(:,1),'r','LineWidth',2.5)
hold on
plot (Shank_Angular_Speed(:,1),'b', 'LineWidth',2.5)
hold on
plot (Thigh_Angular_Speed(:,1), 'k','LineWidth',2.5)
hold on
plot (Trunk_Angular_Speed(:,1), 'g','LineWidth',2.5)
legend ('Foot', 'Shank', 'Thigh','Trunk')


figure('Name','Phase Portrait Normalized')
title('Phase Portrait');
xlabel('Angular Position (º)');
ylabel ('Angular Speed (º/s)');
hold on
plot(Shank_Angle_Normalized(1:TO), Shank_Angular_Speed_Normalized(1:TO), 'g','LineWidth',2.5)
hold on
plot(Shank_Angle_Normalized(TO:end), Shank_Angular_Speed_Normalized(TO:end), 'm','LineWidth',2.5)
hold on
plot(Shank_Angle_Normalized(1),Shank_Angular_Speed_Normalized(1), 'k*', 'LineWidth', 3.0)
hold on
plot(Thigh_Angle_Normalized(1:TO), Thigh_Angular_Speed_Normalized(1:TO), 'b','LineWidth',2.5)
hold on
plot(Thigh_Angle_Normalized(TO:end), Thigh_Angular_Speed_Normalized(TO:end), 'r','LineWidth',2.5)
hold on
plot(Thigh_Angle_Normalized(1),Thigh_Angular_Speed_Normalized(1), 'k*', 'LineWidth', 3.0)
legend ('Shank Contact','Shank Swing', 'TD', 'Thigh Contact', 'Thigh Swing', 'TD')


figure('Name','Phase Portrait')
title('Phase Portrait');
xlabel('Angular Position (º)');
ylabel ('Angular Speed (º/s)');
hold on
plot(Shank_Angle(1:TO), Shank_Angular_Speed(1:TO), 'g','LineWidth',2.5)
hold on
plot(Shank_Angle(TO:end), Shank_Angular_Speed(TO:end), 'm','LineWidth',2.5)
hold on
plot(Shank_Angle(1),Shank_Angular_Speed(1), 'k*', 'LineWidth', 3.0)
hold on
plot(Thigh_Angle(1:TO), Thigh_Angular_Speed(1:TO), 'b','LineWidth',2.5)
hold on
plot(Thigh_Angle(TO:end), Thigh_Angular_Speed(TO:end), 'r','LineWidth',2.5)
hold on
plot(Thigh_Angle(1),Thigh_Angular_Speed(1), 'k*', 'LineWidth', 3.0)
legend ('Shank Contact','Shank Swing', 'TD', 'Thigh Contact','Thigh Swing', 'TD')


figure('Name','Phase Angle and Relative Phase')
title('Phase Angle');
xlabel('Frame (n)');
ylabel ('Phase Angle (º)');
hold on
plot(Shank_Phase_Angle(:,1), 'g','LineWidth',2.5)
hold on
plot(Thigh_Phase_Angle(:,1), 'r','LineWidth',2.5)
hold on
plot (Relative_Phase_Shank_Thigh(:,1), 'b', 'LineWidth', 2.5)
legend('Shank', 'Thigh', 'Relative Phase Shank/Thigh')



figure('Name','Joints Angles')
title('Joints Angles');
hold on
plot(Ankle_Joint_Angle(:,1),'r','LineWidth',2.5)
hold on
plot (Knee_Joint_Angle(:,1),'b', 'LineWidth',2.5)
hold on
plot (Hip_Joint_Angle(:,1), 'k','LineWidth',2.5)
hold on
legend ('Ankle', 'Knee', 'Hip')

figure('Name','Joints Angular Speed')
title('Joints Angles');
hold on
plot(Ankle_Joint_Angular_Speed(:,1),'r','LineWidth',2.5)
hold on
plot (Knee_Joint_Angular_Speed(:,1),'b', 'LineWidth',2.5)
hold on
plot (Hip_Joint_Angular_Speed(:,1), 'k','LineWidth',2.5)
hold on
legend ('Ankle', 'Knee', 'Hip')

%%  SECTION 11
% GROUND REACTION FORCES ESTIMATION %

% Fv = 566.682 - 922.059*(Immersion Ratio) + 3.321 * (Mass) + 176.369 *
% Speed

% Fap = -164.160 + 243.626*(Speed) + 2.443*(Thigh Circumference)
%% SECTION 12 
% ALGORITHM FOR AUTOMATIC TOUCH DOWN AND TAKE OFF EVENTS %
%(O'Connor et al., 2007)

for i = 2:Frames-1
    Velocity_Foot_Algorithm(i) = (Heel(i+1) - Toe(i-1))/(2*dt);
end


%figure ('Name', 'Foot Alghoritm')
%title ('Foot Alghoritm');
%xlabel('Frames (n)');
%ylabel ('Speed (m/s)');
%hold on
%plot(Velocity_Foot_Algorithm,'r');
%hold on
%plot([TD1 TD1],[-0.5 0.5],'-k', 'LineWidth',2.5)
%hold on
%plot([TO TO],[-0.5 0.5],'-g', 'LineWidth',2.5)
%hold on
%plot([TD2-TD1 TD2-TD1],[-0.5 0.5],'-K', 'LineWidth',2.5)

%% SECTION 13 
% EXPORT DATA %

      % SPATIOTEMPORAL DATA        
Export_File_Path_Full_Spatiotemporal = [Export_File_Path '\' File_Name '_SPATIOTEMPORAL' '.xls'];
      
Export_Data_SpatioTemporal = [Speed_Mean Stride_Frequency Stride_Lenght Stride_Duration Contact_Phase_Duration Duty_Factor];
Export_Header_Spatiotemporal = {'Speed_Mean (m/s)', 'Stride_Frequency (Hz)', 'Stride_Lenght (m)', 'Stride_Duration (s)', 'Contact_Phase_Duration (s)', 'Duty Factor (%)'};

Export_Excel_Spatiotemporal = [Export_Header_Spatiotemporal; num2cell(Export_Data_SpatioTemporal)];
xlswrite (Export_File_Path_Full_Spatiotemporal, Export_Excel_Spatiotemporal); 

    
      % ANGULAR DATA
Export_File_Path_Full_Angular = [Export_File_Path '\' File_Name '_ANGULAR' '.xls'];
      
Export_Data_Angular = [Foot_Angle_ROM Shank_Angle_ROM Thigh_Angle_ROM Trunk_Angle_ROM Foot_Angular_Speed_Max Shank_Angular_Speed_Max Thigh_Angular_Speed_Max Trunk_Angular_Speed_Max Ankle_Joint_Angle_ROM Knee_Joint_Angle_ROM Hip_Joint_Angle_ROM Ankle_Joint_Angular_Speed_Max Knee_Joint_Angular_Speed_Max Hip_Joint_Angular_Speed_Max];
Export_Header_Angular = {'Foot_Angle_ROM', 'Shank_Angle_ROM', 'Thigh_Angle_ROM', 'Trunk_Angle_ROM', 'Foot_Angular_Speed_Max', 'Shank_Angular_Speed_Max', 'Thigh_Angular_Speed_Max', 'Trunk_Angular_Speed_Max', 'Ankle_Joint_Angle_ROM', 'Knee_Joint_Angle_ROM', 'Hip_Joint_Angle_ROM', 'Ankle_Joint_Angular_Speed_Max', 'Knee_Joint_Angular_Speed_Max', 'Hip_Joint_Angular_Speed_Max'};

Export_Excel_Angular = [Export_Header_Angular; num2cell(Export_Data_Angular)];
xlswrite (Export_File_Path_Full_Angular, Export_Excel_Angular); 



    % DRAG FORCES DATA
Output_File_Path_Full = [Export_File_Path '\' File_Name '.xls'];
      
Data_DragRMS_SpatioTemporal = [RMS_F_Xiph_Contact RMS_F_Xiph_Swing RMS_F_Umbi_Contact RMS_F_Umbi_Swing RMS_F_UL_Contact RMS_F_UL_Swing RMS_F_LL_Contact RMS_F_LL_Swing RMS_F_FT_Contact RMS_F_FT_Swing Speed_Mean Stride_Frequency Stride_Lenght];
Header_DragRMS_Spatiotemporal = {'RMS_F_Xiph_Contact','RMS_F_Xiph_Swing', 'RMS_F_Umbi_Contact', 'RMS_F_Umbi_Swing' , 'RMS_F_UL_Contact','RMS_F_UL_Swing','RMS_F_LL_Contact', 'RMS_F_LL_Swing', 'RMS_F_FT_Contact', 'RMS_F_FT_Swing', 'Speed_Mean (m/s)', 'Step_Frequency (Hz)', 'Step_Lenght (m)'};

Export_Excel = [Header_DragRMS_Spatiotemporal; num2cell(Data_DragRMS_SpatioTemporal )];
xlswrite (Output_File_Path_Full, Export_Excel); 

%UL_Drag_Curves = [F_UL_Filtered];
%xlswrite('C:\Users\andre\Documents\Mestrado\Congressos\Fisiomecânica 2019\Trabalho Arrasto Parkinson\Material de Apoio\Drag_UL_Data.xls', UL_Drag_Curves);

%LL_Drag_Curves = [F_LL_Filtered];
%xlswrite('C:\Users\andre\Documents\Mestrado\Congressos\Fisiomecânica 2019\Trabalho Arrasto Parkinson\Material de Apoio\Drag_LL_Data.xls', LL_Drag_Curves);

%FT_Drag_Curves = [F_FT_Filtered];
%xlswrite('C:\Users\andre\Documents\Mestrado\Congressos\Fisiomecânica 2019\Trabalho Arrasto Parkinson\Material de Apoio\Drag_FT_Data.xls', FT_Drag_Curves);


  

%% SECTION 14
% GRAPHICS AND GIFS %

[nlTimeStride, ~]=size(Time_Stride_Percentage);
%Time_Stride_Percentage;
F_UL_Filtered_Resized=interpft(F_UL_Filtered, nlTimeStride);
F_LL_Filtered_Resized=interpft(F_LL_Filtered, nlTimeStride);
F_FT_Filtered_Resized=interpft(F_FT_Filtered, nlTimeStride);


figure('Name', 'Drag Force Upper Leg', 'units', 'normalized' ,'outerposition', [0 0 1 1])
title ('Drag Force on Upper Leg during One Stride','FontSize', 48);
xlabel('Stride Percentage (%)','FontSize', 36);
ylabel ('Drag Force (N)','FontSize', 36);
set(gca,'FontSize',30);
hold on
plot(Time_Stride_Percentage(1:TO), F_UL_Filtered_Resized(1:TO), 'b-', 'LineWidth', 4);
hold on;
plot(Time_Stride_Percentage((TO):length(Time_Stride_Percentage)), F_UL_Filtered_Resized((TO):length(Time_Stride_Percentage)), 'r-', 'LineWidth', 4);
legend ('Contact Phase', 'Swing Phase','Location','southwest');
hold off


figure('Name', 'Drag Force Lower Leg', 'units', 'normalized' ,'outerposition', [0 0 1 1])
title ('Drag Force on Lower Leg during One Stride','FontSize', 48);
xlabel('Stride Percentage (%)','FontSize', 36);
ylabel ('Drag Force (N)','FontSize', 36);
set(gca,'FontSize',30);
hold on
plot(Time_Stride_Percentage(1:TO), F_LL_Filtered_Resized(1:TO), 'b-', 'LineWidth', 4);
hold on;
plot(Time_Stride_Percentage(TO:length(Time_Stride_Percentage)), F_LL_Filtered_Resized(TO:length(Time_Stride_Percentage)), 'r-', 'LineWidth', 4);
legend ('Contact Phase', 'Swing Phase','Location','southwest');
hold off

figure('Name', 'Drag Force Foot', 'units', 'normalized' ,'outerposition', [0 0 1 1])
title ('Drag Force on Foot during One Stride','FontSize', 48);
xlabel('Stride Percentage (%)','FontSize', 36);
ylabel ('Drag Force (N)','FontSize', 36);
set(gca,'FontSize',30);
hold on
plot(Time_Stride_Percentage(1:TO), F_FT_Filtered_Resized(1:TO), 'b-', 'LineWidth', 4);
hold on;
plot(Time_Stride_Percentage(TO:length(Time_Stride_Percentage)), F_FT_Filtered_Resized(TO:length(Time_Stride_Percentage)), 'r-', 'LineWidth', 4);
legend ('Contact Phase', 'Swing Phase','Location','southwest');
hold off


%Open Video
%vidObj = VideoReader('C:\Users\andre\Documents\Mestrado\Projetos Pesquisa\Projeto Marcha Água\Piloto\Piloto 06.12.19\Vídeos\Waist\Waist_0.6Step.avi');
%numFrames = get(vidObj,'NumberOfFrames');

%filename = 'DragForceAnimated.gif';
%figure  ('units', 'normalized' ,'outerposition', [0 0 1 1]);
%hold on;
%for k = 1:(length(F_FT_Filtered_Resized)-1)
 %   b = read(vidObj, k);
    
  %  subplot(4,2,[1 2])
   % if k < (TO_in_Percentage+3)
    %    plot(k,F_FT_Filtered_Resized(k), '*b', 'LineWidth', 4), axis([0 100 -10 5]), title('Drag Force on Foot','FontSize', 18), xlabel('Stride Percentage (%)','FontSize', 14), ylabel ('Drag Force (N)','FontSize', 14)
    %    hold on      %hold on function to keep the previous values in the graph
   % else
    %    plot(k,F_FT_Filtered_Resized(k), '*r', 'LineWidth', 4), axis([0 100 -10 5]), title('Drag Force on Foot','FontSize', 18), xlabel('Stride Percentage (%)','FontSize', 14), ylabel ('Drag Force (N)','FontSize', 14)
     %   hold on
   % end
    %hold on
    %subplot(4,2,[5 6 7 8]), imshow(b);
    %drawnow;
    %pause(0000.1);
    
    % gif utilities
    %set(gcf,'color','w'); % set figure background to white
    %drawnow;
    %frame = getframe(1);
    %im = frame2im(frame);
    %[imind,cm] = rgb2ind(im,256);
    %outfile = 'sinewave.gif';
    %if k == 1
     %   imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
   % else
   %     imwrite(imind,cm,filename,'gif','WriteMode','append');
   % end
%end
%hold off