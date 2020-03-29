
% By ANDRÉ IVANISKI MELLO (andreivaniskimello@gmail.com)
% Last edit: 29.03.20


% ROUTINE TO CALCULATE SPATIOTEMPORAL, ANGULAR, COORDINATION, AND KINETICS
% (DRAG FORCE) E GROUND REACTION FORCES FROM KINEMATIC SHALLOW WATER
% WALKING DATA


% INPUTS:
    % POSITION x TIME Matrix of LOWER LIMBR AND TRUNK POINTS
    % SUBJECT BODY MASS
    % IMMERSION DEPTH
    % Touch-Down and Take-Off Frames of ONE STRIDE
    % LOW-PASS BUTTERWORTH FILTER PARAMETERS
    % SAMPLE FREQUENCY
    % NUMBER OF INTERPOLATION POINTS TO NORMALIZE CURVES TO 0 - 100%
    % STRIDE, CONTACT, SWING PHASES DURATIONS
    
% OUTPUTS:
    % SPATIOTEMPORAL
    % ANGULAR (SEGMENTS AND JOINTS)
    % INTRALIMB COORDINATION PARAMETERS
    % DRAG FORCE
    % GROUND REACTION FORCES PEAKS
     
%% DATA ANALYSIS TO INSERT

% Winter residual analysis

%%
%%Orientations%%

% 1º) Insert SubjectName_Depth__Speed_StrideNumber
% 2º) Insert files path of 'Anthropometric' and 'Linear' datas
% 3º) Insert Output file path
% 4º) Insert TD1, TO, TD2 frames
% 5º) Insert Immersion Depth (m) and Subject Mass (kg)
% 6º) Insert Low-Pass Butterworth Filter specifications
% 7º) Insert a number of points (n*100, n is an integer) for interpolation
% of 0-100% curves
% 8º)Run routine


% THIS ROUTINE HAS 15 SECTIONS

% SECTION 01: INFORMATION FOR ROUTINE AND INPUT DATA IMPORTATION
% SECTION 02: DATA INFORMATION PREPARATION 
% SECTION 03: SPATIOTEMPORAL STRIDE PARAMETERS
% SECTION 04: ANGULAR STRIDE PARAMETERS 
% SECTION 05: ANGULAR COORDINATION
% SECTION 06: DRAG FORCE MODEL DATA PREPARATION
% SECTION 07: XIPHOID TRUNK MODEL
% SECTION 08: UMBILICAL TRUNK MODEL
% SECTION 09: UPPER LEG (UL) MODEL
% SECTION 10: LOWER LEG (LL) MODEL
% SECTION 11: FOOT (FT) MODEL
% SECTION 12: DRAG FORCES SUMMATION
% SECTION 13: GROUND REACTION FORCES 
% SECTION 14: ESTIMATION ALGORITHM FOR AUTOMATIC TOUCH DOWN AND TAKE OFF EVENTS 
% SECTION 15: EXPORT DATA GRAPHICS AND GIFS
% SECTION 16: GRAPHICS AND GIFS

clear all
close all
clc

%% SECTION 01
% INFORMATION FOR ROUTINE AND INPUT DATA IMPORTATION %

File_Name =['Lorenzo_Xifoide_0.2_P1'];        %insert 'SubjectName_Depth_Speed_StrideNumber' here, in order to export the final data with the Subject Name_Depth 

%Open Anthropometric Values
Anthropometric_file = fopen('C:\Users\andre\Documents\Mestrado\Projetos Pesquisa\Projeto Marcha Água\Coletas\Data\Cinemetria\Antropometria\Lorenzo_Antropometria.txt');
data = textscan(Anthropometric_file,'%f%s');
fclose(Anthropometric_file);
Anthropometric = data{1,1};

%Open Walk Linear Data
Walk_Linear_file = ('C:\Users\andre\Documents\Mestrado\Projetos Pesquisa\Projeto Marcha Água\Coletas\Data\Cinemetria\Data Out\Lorenzo\Xifoide\Lorenzo_Xifoide_0.2_P1_Positions.txt');
delimiterIn = '\t';
Linear_Data = importdata(Walk_Linear_file,delimiterIn);
Linear_Data = Linear_Data.';

%Insert Export (Output) File Path here
Export_File_Path = ['C:\Users\andre\Documents\Mestrado\Projetos Pesquisa\Projeto Marcha Água\Coletas\Data\Cinemetria\Data Out\Lorenzo\Xifoide\Planilhas_MatLab'];     

% Insert TD and TO frames 
    % TD: First frame that the foot touch the ground
    % TO: First frame that the foot is out of the ground

TD1 = 39;
TO = 146;                              
TD2 = 247;

TO = TO - TD1;     % This is adjusting TO frame number considering TD1 as Frame 1
TD2 = TD2 - TD1;   % This is adjusting TD2 frame number considering TD1 as Frame 1 


% Immersion depth (m)
Immersion_depth = 1.2;


% Subject mass (kg)
Mass = 75 ;


% Low-Pass Butterworth Filter Specifications
fsample = 60;
dt=1/fsample;
fcut=4;
order=2;

% Information for Curves Interpolation and Resizing to 0 - 100% of Stride
% Duration

Interpolation_Points = 300;                         % Set a 100*n value, higher than the total Stride Frames number! n is an Integer.

%% SECTION 02
% DATA INFORMATION PREPARATION %

    % Body Weight
Body_Weight = Mass * 9.8;


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


Linear_Data = Linear_Data(TD1:(TD1+TD2-1),:);

Frames = size(Linear_Data,1);                    % Calculate the total frames number
Toe_Raw = Linear_Data(:,1:2);                     
Heel_Raw = Linear_Data(:,3:4);
Ankle_Raw = Linear_Data(:,5:6);
Knee_Raw = Linear_Data(:,7:8);
Trochanter_Raw = Linear_Data(:,9:10);
Trunk_Point_Raw = Linear_Data(:,11:12);

Frames_Vector = [1:1:Frames];
Frames_Vector = Frames_Vector.';
Vector_Interpolation_Stride = linspace(1, Frames, Interpolation_Points);
Vector_Interpolation_Contact = linspace(1, (TO-1), Interpolation_Points);
Vector_Interpolation_Swing = linspace(1, (TD2-TO), Interpolation_Points);


%%% LOWPASS BUTTERWORTH FILTER %%%

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
Time_Vector = [0:dt:(TD2/fsample)];                              %Create a Time Column (Vector)
Time_Stride_Percentage = Time_Vector/(max(Time_Vector))*100;     %Create a Time Vector (row) (1xn) in Percentage of the Stride
Time_Stride_Percentage = Time_Stride_Percentage.';               %Rotate Time Vector (column) in Percentage of the Stride to (nx1)
[nlTimeStride, ~]=size(Time_Stride_Percentage);
TO_in_Percentage_Decimal = TO/(nlTimeStride-1)*100;
TO_in_Percentage = fix(TO_in_Percentage_Decimal);                 %Calculate TO in Percentage of the Stride

Stride_Lenght = (Heel_x(end) - Heel_x(1));                        %Final position - Initial Position
Speed_Mean = (Stride_Lenght/Time_Total);                          %Speed_Mean: is the mean speed value over the entire walk.
Stride_Frequency = Speed_Mean/Stride_Lenght;

Stride_Duration = (TD2)*dt;
Contact_Phase_Duration = (TO)*dt;
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


            % INTERPOLATION TO 0-100% OF STRIDE DURATION (FULL STRIDE,
            % CONTACT AND SWING PHASES)
           
% Foot
Foot_Angle_Interpolation_Stride = interp1(Foot_Angle, Vector_Interpolation_Stride, 'linear');
Foot_Angle_Interpolation_Stride = Foot_Angle_Interpolation_Stride';
Foot_Angle_Interpolation_Stride = reshape (Foot_Angle_Interpolation_Stride, [], 100);
Foot_Angle_100_Stride = mean (Foot_Angle_Interpolation_Stride);
Foot_Angle_100_Stride = Foot_Angle_100_Stride.';

Foot_Angle_Interpolation_Contact = interp1(Foot_Angle(1:TO-1), Vector_Interpolation_Contact, 'linear');
Foot_Angle_Interpolation_Contact = Foot_Angle_Interpolation_Contact';
Foot_Angle_Interpolation_Contact = reshape (Foot_Angle_Interpolation_Contact, [], 100);
Foot_Angle_100_Contact = mean (Foot_Angle_Interpolation_Contact);
Foot_Angle_100_Contact = Foot_Angle_100_Contact.';

Foot_Angle_Interpolation_Swing = interp1(Foot_Angle(TO:end), Vector_Interpolation_Swing, 'linear');
Foot_Angle_Interpolation_Swing = Foot_Angle_Interpolation_Swing';
Foot_Angle_Interpolation_Swing = reshape (Foot_Angle_Interpolation_Swing, [], 100);
Foot_Angle_100_Swing = mean (Foot_Angle_Interpolation_Swing);
Foot_Angle_100_Swing = Foot_Angle_100_Swing.';

% Shank
Shank_Angle_Interpolation_Stride = interp1(Shank_Angle, Vector_Interpolation_Stride, 'linear');
Shank_Angle_Interpolation_Stride = Shank_Angle_Interpolation_Stride';
Shank_Angle_Interpolation_Stride = reshape (Shank_Angle_Interpolation_Stride, [], 100);
Shank_Angle_100_Stride = mean (Shank_Angle_Interpolation_Stride);
Shank_Angle_100_Stride = Shank_Angle_100_Stride.';

Shank_Angle_Interpolation_Contact = interp1(Shank_Angle(1:TO-1), Vector_Interpolation_Contact, 'linear');
Shank_Angle_Interpolation_Contact = Shank_Angle_Interpolation_Contact';
Shank_Angle_Interpolation_Contact = reshape (Shank_Angle_Interpolation_Contact, [], 100);
Shank_Angle_100_Contact = mean (Shank_Angle_Interpolation_Contact);
Shank_Angle_100_Contact = Shank_Angle_100_Contact.';

Shank_Angle_Interpolation_Swing = interp1(Shank_Angle(TO:end), Vector_Interpolation_Swing, 'linear');
Shank_Angle_Interpolation_Swing = Shank_Angle_Interpolation_Swing';
Shank_Angle_Interpolation_Swing = reshape (Shank_Angle_Interpolation_Swing, [], 100);
Shank_Angle_100_Swing = mean (Shank_Angle_Interpolation_Swing);
Shank_Angle_100_Swing = Shank_Angle_100_Swing.';

% Thigh
Thigh_Angle_Interpolation_Stride = interp1(Thigh_Angle, Vector_Interpolation_Stride, 'linear');
Thigh_Angle_Interpolation_Stride = Thigh_Angle_Interpolation_Stride';
Thigh_Angle_Interpolation_Stride = reshape (Thigh_Angle_Interpolation_Stride, [], 100);
Thigh_Angle_100_Stride = mean (Thigh_Angle_Interpolation_Stride);
Thigh_Angle_100_Stride = Thigh_Angle_100_Stride.';

Thigh_Angle_Interpolation_Contact = interp1(Thigh_Angle(1:TO-1), Vector_Interpolation_Contact, 'linear');
Thigh_Angle_Interpolation_Contact = Thigh_Angle_Interpolation_Contact';
Thigh_Angle_Interpolation_Contact = reshape (Thigh_Angle_Interpolation_Contact, [], 100);
Thigh_Angle_100_Contact = mean (Thigh_Angle_Interpolation_Contact);
Thigh_Angle_100_Contact = Thigh_Angle_100_Contact.';

Thigh_Angle_Interpolation_Swing = interp1(Thigh_Angle(TO:end), Vector_Interpolation_Swing, 'linear');
Thigh_Angle_Interpolation_Swing = Thigh_Angle_Interpolation_Swing';
Thigh_Angle_Interpolation_Swing = reshape (Thigh_Angle_Interpolation_Swing, [], 100);
Thigh_Angle_100_Swing = mean (Thigh_Angle_Interpolation_Swing);
Thigh_Angle_100_Swing = Thigh_Angle_100_Swing.';

% Trunk
Trunk_Angle_Interpolation_Stride = interp1(Trunk_Angle, Vector_Interpolation_Stride, 'linear');
Trunk_Angle_Interpolation_Stride = Trunk_Angle_Interpolation_Stride';
Trunk_Angle_Interpolation_Stride = reshape (Trunk_Angle_Interpolation_Stride, [], 100);
Trunk_Angle_100_Stride = mean (Trunk_Angle_Interpolation_Stride);
Trunk_Angle_100_Stride = Trunk_Angle_100_Stride.';

Trunk_Angle_Interpolation_Contact = interp1(Trunk_Angle(1:TO-1), Vector_Interpolation_Contact, 'linear');
Trunk_Angle_Interpolation_Contact = Trunk_Angle_Interpolation_Contact';
Trunk_Angle_Interpolation_Contact = reshape (Trunk_Angle_Interpolation_Contact, [], 100);
Trunk_Angle_100_Contact = mean (Trunk_Angle_Interpolation_Contact);
Trunk_Angle_100_Contact = Trunk_Angle_100_Contact.';

Trunk_Angle_Interpolation_Swing = interp1(Trunk_Angle(TO:end), Vector_Interpolation_Swing, 'linear');
Trunk_Angle_Interpolation_Swing = Trunk_Angle_Interpolation_Swing';
Trunk_Angle_Interpolation_Swing = reshape (Trunk_Angle_Interpolation_Swing, [], 100);
Trunk_Angle_100_Swing = mean (Trunk_Angle_Interpolation_Swing);
Trunk_Angle_100_Swing = Trunk_Angle_100_Swing.';


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


            % INTERPOLATION TO 0-100% OF STRIDE DURATION (FULL STRIDE,
            % CONTACT AND SWING PHASES)

% Foot
for i = 2:100-1
    Foot_Angular_Speed_100_Stride (i) = (Foot_Angle_100_Stride(i+1) - Foot_Angle_100_Stride(i-1))/(2*dt);
end
Foot_Angular_Speed_100_Stride = interpft(Foot_Angular_Speed_100_Stride,100);
Foot_Angular_Speed_100_Stride = Foot_Angular_Speed_100_Stride';
Foot_Angular_Speed_100_Stride_Max = max(Foot_Angular_Speed_100_Stride);

% Shank
for i = 2:100-1
    Shank_Angular_Speed_100_Stride (i) = (Shank_Angle_100_Stride(i+1) - Shank_Angle_100_Stride(i-1))/(2*dt);
end
Shank_Angular_Speed_100_Stride = interpft(Shank_Angular_Speed_100_Stride,100);
Shank_Angular_Speed_100_Stride = Shank_Angular_Speed_100_Stride';
Shank_Angular_Speed_100_Stride_Max = max(Shank_Angular_Speed_100_Stride);

% Thigh
for i = 2:100-1
    Thigh_Angular_Speed_100_Stride (i) = (Thigh_Angle_100_Stride(i+1) - Thigh_Angle_100_Stride(i-1))/(2*dt);
end
Thigh_Angular_Speed_100_Stride = interpft(Thigh_Angular_Speed_100_Stride,100);
Thigh_Angular_Speed_100_Stride = Thigh_Angular_Speed_100_Stride';
Thigh_Angular_Speed_100_Stride_Max = max(Thigh_Angular_Speed_100_Stride);



    % JOINT'S ANGULAR POSITION

Ankle_Joint_Angle = Shank_Angle - Foot_Angle;
Knee_Joint_Angle = Thigh_Angle - Shank_Angle;
Hip_Joint_Angle = Trunk_Angle - Thigh_Angle;

    
             % INTERPOLATION TO 0-100% OF STRIDE DURATION (FULL STRIDE,
            % CONTACT AND SWING PHASES)
            
% Ankle_Joint
Ankle_Joint_Angle_Interpolation_Stride = interp1(Ankle_Joint_Angle, Vector_Interpolation_Stride, 'linear');
Ankle_Joint_Angle_Interpolation_Stride = Ankle_Joint_Angle_Interpolation_Stride';
Ankle_Joint_Angle_Interpolation_Stride = reshape (Ankle_Joint_Angle_Interpolation_Stride, [], 100);
Ankle_Joint_Angle_100_Stride = mean (Ankle_Joint_Angle_Interpolation_Stride);
Ankle_Joint_Angle_100_Stride = Ankle_Joint_Angle_100_Stride.';

Ankle_Joint_Angle_Interpolation_Contact = interp1(Ankle_Joint_Angle(1:TO-1), Vector_Interpolation_Contact, 'linear');
Ankle_Joint_Angle_Interpolation_Contact = Ankle_Joint_Angle_Interpolation_Contact';
Ankle_Joint_Angle_Interpolation_Contact = reshape (Ankle_Joint_Angle_Interpolation_Contact, [], 100);
Ankle_Joint_Angle_100_Contact = mean (Ankle_Joint_Angle_Interpolation_Contact);
Ankle_Joint_Angle_100_Contact = Ankle_Joint_Angle_100_Contact.';

Ankle_Joint_Angle_Interpolation_Swing = interp1(Ankle_Joint_Angle(TO:end), Vector_Interpolation_Swing, 'linear');
Ankle_Joint_Angle_Interpolation_Swing = Ankle_Joint_Angle_Interpolation_Swing';
Ankle_Joint_Angle_Interpolation_Swing = reshape (Ankle_Joint_Angle_Interpolation_Swing, [], 100);
Ankle_Joint_Angle_100_Swing = mean (Ankle_Joint_Angle_Interpolation_Swing);
Ankle_Joint_Angle_100_Swing = Ankle_Joint_Angle_100_Swing.';

% Knee_Joint
Knee_Joint_Angle_Interpolation_Stride = interp1(Knee_Joint_Angle, Vector_Interpolation_Stride, 'linear');
Knee_Joint_Angle_Interpolation_Stride = Knee_Joint_Angle_Interpolation_Stride';
Knee_Joint_Angle_Interpolation_Stride = reshape (Knee_Joint_Angle_Interpolation_Stride, [], 100);
Knee_Joint_Angle_100_Stride = mean (Knee_Joint_Angle_Interpolation_Stride);
Knee_Joint_Angle_100_Stride = Knee_Joint_Angle_100_Stride.';

Knee_Joint_Angle_Interpolation_Contact = interp1(Knee_Joint_Angle(1:TO-1), Vector_Interpolation_Contact, 'linear');
Knee_Joint_Angle_Interpolation_Contact = Knee_Joint_Angle_Interpolation_Contact';
Knee_Joint_Angle_Interpolation_Contact = reshape (Knee_Joint_Angle_Interpolation_Contact, [], 100);
Knee_Joint_Angle_100_Contact = mean (Knee_Joint_Angle_Interpolation_Contact);
Knee_Joint_Angle_100_Contact = Knee_Joint_Angle_100_Contact.';

Knee_Joint_Angle_Interpolation_Swing = interp1(Knee_Joint_Angle(TO:end), Vector_Interpolation_Swing, 'linear');
Knee_Joint_Angle_Interpolation_Swing = Knee_Joint_Angle_Interpolation_Swing';
Knee_Joint_Angle_Interpolation_Swing = reshape (Knee_Joint_Angle_Interpolation_Swing, [], 100);
Knee_Joint_Angle_100_Swing = mean (Knee_Joint_Angle_Interpolation_Swing);
Knee_Joint_Angle_100_Swing = Knee_Joint_Angle_100_Swing.';

% Hip_Joint
Hip_Joint_Angle_Interpolation_Stride = interp1(Hip_Joint_Angle, Vector_Interpolation_Stride, 'linear');
Hip_Joint_Angle_Interpolation_Stride = Hip_Joint_Angle_Interpolation_Stride';
Hip_Joint_Angle_Interpolation_Stride = reshape (Hip_Joint_Angle_Interpolation_Stride, [], 100);
Hip_Joint_Angle_100_Stride = mean (Hip_Joint_Angle_Interpolation_Stride);
Hip_Joint_Angle_100_Stride = Hip_Joint_Angle_100_Stride.';

Hip_Joint_Angle_Interpolation_Contact = interp1(Hip_Joint_Angle(1:TO-1), Vector_Interpolation_Contact, 'linear');
Hip_Joint_Angle_Interpolation_Contact = Hip_Joint_Angle_Interpolation_Contact';
Hip_Joint_Angle_Interpolation_Contact = reshape (Hip_Joint_Angle_Interpolation_Contact, [], 100);
Hip_Joint_Angle_100_Contact = mean (Hip_Joint_Angle_Interpolation_Contact);
Hip_Joint_Angle_100_Contact = Hip_Joint_Angle_100_Contact.';

Hip_Joint_Angle_Interpolation_Swing = interp1(Hip_Joint_Angle(TO:end), Vector_Interpolation_Swing, 'linear');
Hip_Joint_Angle_Interpolation_Swing = Hip_Joint_Angle_Interpolation_Swing';
Hip_Joint_Angle_Interpolation_Swing = reshape (Hip_Joint_Angle_Interpolation_Swing, [], 100);
Hip_Joint_Angle_100_Swing = mean (Hip_Joint_Angle_Interpolation_Swing);
Hip_Joint_Angle_100_Swing = Hip_Joint_Angle_100_Swing.';


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

%% SECTION 05
% ANGULAR COORDINATION %

    % PHASE PORTRAIT
        % (Degani e Danna-dos-Santo, 2007)
       
% Foot
Foot_Angle_100_Stride_Normalized =  Foot_Angle_100_Stride/max(Foot_Angle_100_Stride(:,1));
Foot_Angular_Speed_100_Stride_Normalized  = Foot_Angular_Speed_100_Stride/max(Foot_Angular_Speed_100_Stride(:,1));
Foot_Phase_Angle = atand(Foot_Angular_Speed_100_Stride_Normalized(:,1)./Foot_Angle_100_Stride_Normalized (:,1));      

% Shank
Shank_Angle_100_Stride_Normalized = Shank_Angle_100_Stride/max(Shank_Angle_100_Stride(:,1));
Shank_Angular_Speed_100_Stride_Normalized  = Shank_Angular_Speed_100_Stride/max(Shank_Angular_Speed_100_Stride(:,1));
Shank_Phase_Angle = atand(Shank_Angular_Speed_100_Stride_Normalized(:,1)./Shank_Angle_100_Stride_Normalized(:,1)); 

% Thigh
Thigh_Angle_100_Stride_Normalized = Thigh_Angle_100_Stride/max(Thigh_Angle_100_Stride(:,1));
Thigh_Angular_Speed_100_Stride_Normalized  = Thigh_Angular_Speed_100_Stride/max(Thigh_Angular_Speed_100_Stride(:,1));
Thigh_Phase_Angle = atand(Thigh_Angular_Speed_100_Stride_Normalized(:,1)./Thigh_Angle_100_Stride_Normalized(:,1)); 


    % RELATIVE PHASE
Relative_Phase_Foot_Shank =  Foot_Phase_Angle - Shank_Phase_Angle;
Relative_Phase_Shank_Thigh = Shank_Phase_Angle - Thigh_Phase_Angle;


    % VECTOR CODING ANALYSIS
        % (Chang et al. 2008; Silvernail et al. 2018; Celestino et al. 2019)

Foot_Angle_100_Stride_diff = diff(Foot_Angle_100_Stride);
Foot_Angle_100_Contact_diff =diff(Foot_Angle_100_Contact);
Foot_Angle_100_Swing_diff = diff(Foot_Angle_100_Swing);
        

Shank_Angle_100_Stride_diff = diff(Shank_Angle_100_Stride);
Shank_Angle_100_Contact_diff =diff(Shank_Angle_100_Contact);
Shank_Angle_100_Swing_diff = diff(Shank_Angle_100_Swing);
    
Thigh_Angle_100_Stride_diff = diff(Thigh_Angle_100_Stride);
Thigh_Angle_100_Contact_diff =diff(Thigh_Angle_100_Contact);
Thigh_Angle_100_Swing_diff = diff(Thigh_Angle_100_Swing);

Coupling_Angle_Shank_Foot_Stride = interpft(atand(Foot_Angle_100_Stride_diff./Shank_Angle_100_Stride_diff),100);
Coupling_Angle_Shank_Foot_Contact = interpft(atand(Foot_Angle_100_Contact_diff./Shank_Angle_100_Contact_diff),100);     
Coupling_Angle_Shank_Foot_Swing = interpft(atand(Foot_Angle_100_Swing_diff./Shank_Angle_100_Swing_diff),100);

Coupling_Angle_Thigh_Shank_Stride = interpft(atand(Shank_Angle_100_Stride_diff./Thigh_Angle_100_Stride_diff),100);
Coupling_Angle_Thigh_Shank_Contact = interpft(atand(Shank_Angle_100_Contact_diff./Thigh_Angle_100_Contact_diff),100);
Coupling_Angle_Thigh_Shank_Swing = interpft(atand(Shank_Angle_100_Swing_diff./Thigh_Angle_100_Swing_diff),100);

%% SECTION 06
% DRAG FORCE MODEL DATA PREPARATION %
% (Orselli & Duarte, 2011)

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

%% GENERAL FUNCTIONS FOR DRAG FORCE %%%

%%dA(z) = 2*((Ad-Ap)z + Ap*L)/L))*dz
%%v(z)=((z*Vd)+(L-z)*Vp)/L
%%dAP(z) = d(z).sen(angle)   , dAP: Projected Area at Vertical Axis
%%(perpendicular to Velocity Vector)

%% SECTION 07
% XIPHOID TRUNK MODEL %

%Ap = HipBreadth/2
%Ad = XiphoidC/2(pi)
%L = SittingH - (Stature - SubsternalH)

% XIPHOID TRUNK SIZE MODELING (CYLINDER)
Ap_Xiph = HipB/2;
Ad_Xiph = XiphoidC/(2*pi);
L_Xiph = SittingH - (Stature - SubsternalH);


% DIVISION OF SEGMENT AREA IN dZ PARTS

% dA(z) = 2*((Ad-Ap)z + Ap*L)/L))*dz
% z = (0+dL_Xiph(i))

%Xiphoid Trunk length differential
dL_Xiph = 0:0.001:max(L_Xiph);
 

%dA_Xiph = dA(z)
for i=1:length(dL_Xiph)
    dA_Xiph(i) =(((Ad_Xiph-Ap_Xiph)*(0+dL_Xiph(i))+(Ap_Xiph*L_Xiph))/L_Xiph)*0.001*2;
end


% TOTAL SEGMENT FRONTAL AREA
A_Xiph = sum(dA_Xiph);


% FRONTAL PROJECTED AREA OF EACH dZ PART AREA

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


% LINEAR VELOCITY OF EACH dZ PART

%%v(z)=((z*Vd)+(L-z)*Vp)/L
%Xiph Velocity
%Vp: Velocity_Throcanter
%Vd: Velocity_Xiphoid
%V_Xiph = v(z)

% NESTED FOR LOOP: 
% FIRST FOR LOOP IS TO TIME VECTOR
% SECOND FOR LOOP IS TO SEGMENT dL(z)

%V_Xiph [frames dL_Xiph]

V_Xiph = zeros(length(Velocity_Throcanter), length(dL_Xiph));   %preallocation: to optimize routine running time

for j=1:1:length(Velocity_Throcanter)
    
    for q=1:length(dL_Xiph)
        V_Xiph (j,q) = (((0+dL_Xiph(q))*Velocity_Trunk(j))+((L_Xiph-(0+dL_Xiph(q)))*Velocity_Throcanter(j)))/L_Xiph;
    end
end


% DRAG FORCE IN EACH dZ PART

%Cd: drag coefficient;
%W: water density (values from Orselli & Duarte, 2011)
Cd = 1;
W = 1000;

for j=1:size(V_Xiph,1)
    for q=1:size(V_Xiph,2)
        dF_Xiph(j,q) = Cd*W*0.5*(V_Xiph(j,q).^2)*(dAP_Xiph(j,q));
    end
end


% TOTAL DRAG FORCE ON SEGMENT: SUM OF DRAG FORCE OF EACH dZ PART
F_Xiph = (-1)*sum(dF_Xiph,2);


% FILTERING OF DRAG FORCE
fcut=8;
order=2;
F_Xiph = matfiltfilt(dt, fcut, order, F_Xiph);


% DIVISION OF DRAG FORCE IN CONTACT AND SWING PHASES
F_Xiph_Contact = F_Xiph(1:TO-1);
F_Xiph_Swing = F_Xiph(TO:length(F_Xiph));

F_Xiph_Contact_Max = max(F_Xiph_Contact);
F_Xiph_Contact_Min = min(F_Xiph_Contact);
F_Xiph_Contact_Mean = mean(F_Xiph_Contact);

F_Xiph_Swing_Max = max(F_Xiph_Swing);
F_Xiph_Swing_Min = min(F_Xiph_Swing);
F_Xiph_Swing_Mean = mean(F_Xiph_Swing);


% CALCULATION OF RMS DRAG FORCE OF CONTACT AND SWING PHASES
RMS_F_Xiph_Contact = rms(F_Xiph_Contact);
RMS_F_Xiph_Swing = rms(F_Xiph_Swing);


% NORMALIZATION OF RMS DRAG FORCE BY SEGMENT FRONTAL AREA
RMS_F_Xiph_AreaNormalized_Contact = RMS_F_Xiph_Contact/A_Xiph; 
RMS_F_Xiph_AreaNormalized_Swing = RMS_F_Xiph_Swing/A_Xiph; 

F_Xiph_AreaNormalized = F_Xiph/A_Xiph;


% NORMALIZATION OF DRAG FORCE BY BODY WEIGHT 
F_Xiph_WeightNormalized = (F_Xiph/Body_Weight)*100;

F_Xiph_WeightNormalized_Contact = F_Xiph_WeightNormalized(1:TO-1);
F_Xiph_WeightNormalized_Swing = F_Xiph_WeightNormalized(TO:length(F_Xiph_WeightNormalized));

RMS_F_Xiph_WeightNormalized_Contact = rms(F_Xiph_WeightNormalized_Contact);
RMS_F_Xiph_WeightNormalized_Swing = rms(F_Xiph_WeightNormalized_Swing);


        % INTERPOLATION TO 0-100% OF STRIDE DURATION
        
% DRAG FORCE XIPHOID
F_Xiph_Interpolation_Stride = interp1(F_Xiph, Vector_Interpolation_Stride, 'linear');
F_Xiph_Interpolation_Stride = F_Xiph_Interpolation_Stride';
F_Xiph_Interpolation_Stride = reshape (F_Xiph_Interpolation_Stride, [], 100);
F_Xiph_100_Stride = mean(F_Xiph_Interpolation_Stride);
F_Xiph_100_Stride = F_Xiph_100_Stride.';

F_Xiph_Interpolation_Contact = interp1(F_Xiph(1:TO-1), Vector_Interpolation_Contact, 'linear');
F_Xiph_Interpolation_Contact = F_Xiph_Interpolation_Contact';
F_Xiph_Interpolation_Contact = reshape (F_Xiph_Interpolation_Contact, [], 100);
F_Xiph_100_Contact = mean (F_Xiph_Interpolation_Contact);
F_Xiph_100_Contact = F_Xiph_100_Contact.';

F_Xiph_Interpolation_Swing = interp1(F_Xiph(TO:end), Vector_Interpolation_Swing, 'linear');
F_Xiph_Interpolation_Swing = F_Xiph_Interpolation_Swing';
F_Xiph_Interpolation_Swing = reshape (F_Xiph_Interpolation_Swing, [], 100);
F_Xiph_100_Swing = mean (F_Xiph_Interpolation_Swing);
F_Xiph_100_Swing = F_Xiph_100_Swing.';


% DRAG FORCE XIPHOID NORMALIZED BY BODY  WEIGHT
F_Xiph_WeightNormalized_Interpolation_Stride = interp1(F_Xiph_WeightNormalized, Vector_Interpolation_Stride, 'linear');
F_Xiph_WeightNormalized_Interpolation_Stride = F_Xiph_WeightNormalized_Interpolation_Stride';
F_Xiph_WeightNormalized_Interpolation_Stride = reshape (F_Xiph_WeightNormalized_Interpolation_Stride, [], 100);
F_Xiph_WeightNormalized_100_Stride = mean(F_Xiph_WeightNormalized_Interpolation_Stride);
F_Xiph_WeightNormalized_100_Stride = F_Xiph_WeightNormalized_100_Stride.';

F_Xiph_WeightNormalized_Interpolation_Contact = interp1(F_Xiph_WeightNormalized(1:TO-1), Vector_Interpolation_Contact, 'linear');
F_Xiph_WeightNormalized_Interpolation_Contact = F_Xiph_WeightNormalized_Interpolation_Contact';
F_Xiph_WeightNormalized_Interpolation_Contact = reshape (F_Xiph_WeightNormalized_Interpolation_Contact, [], 100);
F_Xiph_WeightNormalized_100_Contact = mean (F_Xiph_WeightNormalized_Interpolation_Contact);
F_Xiph_WeightNormalized_100_Contact = F_Xiph_WeightNormalized_100_Contact.';

F_Xiph_WeightNormalized_Interpolation_Swing = interp1(F_Xiph_WeightNormalized(TO:end), Vector_Interpolation_Swing, 'linear');
F_Xiph_WeightNormalized_Interpolation_Swing = F_Xiph_WeightNormalized_Interpolation_Swing';
F_Xiph_WeightNormalized_Interpolation_Swing = reshape (F_Xiph_WeightNormalized_Interpolation_Swing, [], 100);
F_Xiph_WeightNormalized_100_Swing = mean (F_Xiph_WeightNormalized_Interpolation_Swing);
F_Xiph_WeightNormalized_100_Swing = F_Xiph_WeightNormalized_100_Swing.';


% DRAG FORCE XIPHOID NORMALIZED BY SEGMENT FRONTAL AREA
F_Xiph_AreaNormalized_Interpolation_Stride = interp1(F_Xiph_AreaNormalized, Vector_Interpolation_Stride, 'linear');
F_Xiph_AreaNormalized_Interpolation_Stride = F_Xiph_AreaNormalized_Interpolation_Stride';
F_Xiph_AreaNormalized_Interpolation_Stride = reshape (F_Xiph_AreaNormalized_Interpolation_Stride, [], 100);
F_Xiph_AreaNormalized_100_Stride = mean(F_Xiph_AreaNormalized_Interpolation_Stride);
F_Xiph_AreaNormalized_100_Stride = F_Xiph_AreaNormalized_100_Stride.';

F_Xiph_AreaNormalized_Interpolation_Contact = interp1(F_Xiph_AreaNormalized(1:TO-1), Vector_Interpolation_Contact, 'linear');
F_Xiph_AreaNormalized_Interpolation_Contact = F_Xiph_AreaNormalized_Interpolation_Contact';
F_Xiph_AreaNormalized_Interpolation_Contact = reshape (F_Xiph_AreaNormalized_Interpolation_Contact, [], 100);
F_Xiph_AreaNormalized_100_Contact = mean (F_Xiph_AreaNormalized_Interpolation_Contact);
F_Xiph_AreaNormalized_100_Contact = F_Xiph_AreaNormalized_100_Contact.';

F_Xiph_AreaNormalized_Interpolation_Swing = interp1(F_Xiph_AreaNormalized(TO:end), Vector_Interpolation_Swing, 'linear');
F_Xiph_AreaNormalized_Interpolation_Swing = F_Xiph_AreaNormalized_Interpolation_Swing';
F_Xiph_AreaNormalized_Interpolation_Swing = reshape (F_Xiph_AreaNormalized_Interpolation_Swing, [], 100);
F_Xiph_AreaNormalized_100_Swing = mean (F_Xiph_AreaNormalized_Interpolation_Swing);
F_Xiph_AreaNormalized_100_Swing = F_Xiph_AreaNormalized_100_Swing.';

%% SECTION 08
% UMBILICAL TRUNK MODEL %

%Ap = HipBreadth/2
%Ad = UmbilicalC/2(pi)
%L = SittingH - (Stature - UmbilicalH)

% UMBILICAL TRUNK SIZE MODELING (CYLINDER)
Ap_Umbi = HipB/2;
Ad_Umbi = UmbilicalC/(2*pi);
L_Umbi = SittingH - (Stature - UmbilicalH);


% DIVISION OF SEGMENT AREA IN dZ PARTS

% dA(z) = 2*((Ad-Ap)z + Ap*L)/L))*dz
% z = (0+dL_UL(i))

%Umbilical Trunk length differential
dL_Umbi = 0:0.001:max(L_Umbi);

%dA_Umbi = dA(z)
for i=1:length(dL_Umbi)
    dA_Umbi(i) =(((Ad_Umbi-Ap_Umbi)*(0+dL_Umbi(i))+(Ap_Umbi*L_Umbi))/L_Umbi)*0.001*2;
end


% TOTAL SEGMENT FRONTAL AREA
A_Umbi = sum(dA_Umbi);


% FRONTAL PROJECTED AREA OF EACH dZ PART AREA

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


% LINEAR VELOCITY OF EACH dZ PART

%%v(z)=((z*Vd)+(L-z)*Vp)/L
%Umbi Velocity
%Vp: Velocity_Throcanter
%Vd: Velocity_Umbilical
%V_Umbi = v(z)

% NESTED FOR LOOP: 
% FIRST FOR LOOP IS TO TIME VECTOR
% SECOND FOR LOOP IS TO SEGMENT dL(z)

%V_Umbi [frames dL_Umbi]

V_Umbi = zeros(length(Velocity_Throcanter), length(dL_Umbi));   %preallocation: to optimize routine running time

for j=1:1:length(Velocity_Throcanter)
    for q=1:length(dL_Umbi)
        V_Umbi (j,q) = (((0+dL_Umbi(q))*Velocity_Trunk(j))+((L_Umbi-(0+dL_Umbi(q)))*Velocity_Throcanter(j)))/L_Umbi;
    end
end


% DRAG FORCE IN EACH dZ PART

%Cd: drag coefficient;
%W: water density (values from Orselli & Duarte, 2011)
Cd = 1;
W = 1000;

for j=1:size(V_Umbi,1)
    for q=1:size(V_Umbi,2)
        dF_Umbi(j,q) = Cd*W*0.5*(V_Umbi(j,q).^2)*(dAP_Umbi(j,q));
    end
end


% TOTAL DRAG FORCE ON SEGMENT: SUM OF DRAG FORCE OF EACH dZ PART
F_Umbi = (-1)*sum(dF_Umbi,2);


% FILTERING OF DRAG FORCE
fcut=8;
order=2;
F_Umbi = matfiltfilt(dt, fcut, order, F_Umbi);


% DIVISION OF DRAG FORCE IN CONTACT AND SWING PHASES
F_Umbi_Contact = F_Umbi(1:TO-1);
F_Umbi_Swing = F_Umbi(TO:length(F_Umbi));

F_Umbi_Contact_Max = max(F_Umbi_Contact);
F_Umbi_Contact_Min = min(F_Umbi_Contact);
F_Umbi_Contact_Mean = mean(F_Umbi_Contact);

F_Umbi_Swing_Max = max(F_Umbi_Swing);
F_Umbi_Swing_Min = min(F_Umbi_Swing);
F_Umbi_Swing_Mean = mean(F_Umbi_Swing);


% CALCULATION OF RMS DRAG FORCE OF CONTACT AND SWING PHASES
RMS_F_Umbi_Contact = rms(F_Umbi_Contact);
RMS_F_Umbi_Swing = rms(F_Umbi_Swing);


% NORMALIZATION OF RMS DRAG FORCE BY SEGMENT FRONTAL AREA
RMS_F_Umbi_AreaNormalized_Contact = RMS_F_Umbi_Contact/A_Umbi; 
RMS_F_Umbi_AreaNormalized_Swing = RMS_F_Umbi_Swing/A_Umbi; 

F_Umbi_AreaNormalized = F_Umbi/A_Umbi;

% NORMALIZATION OF DRAG FORCE BY BODY WEIGHT 
F_Umbi_WeightNormalized = (F_Umbi/Body_Weight)*100;

F_Umbi_WeightNormalized_Contact = F_Umbi_WeightNormalized(1:TO-1);
F_Umbi_WeightNormalized_Swing = F_Umbi_WeightNormalized(TO:length(F_Umbi_WeightNormalized));

RMS_F_Umbi_WeightNormalized_Contact = rms(F_Umbi_WeightNormalized_Contact);
RMS_F_Umbi_WeightNormalized_Swing = rms(F_Umbi_WeightNormalized_Swing);


        % INTERPOLATION TO 0-100% OF STRIDE DURATION
        
% DRAG FORCE UMBILICAL
F_Umbi_Interpolation_Stride = interp1(F_Umbi, Vector_Interpolation_Stride, 'linear');
F_Umbi_Interpolation_Stride = F_Umbi_Interpolation_Stride';
F_Umbi_Interpolation_Stride = reshape (F_Umbi_Interpolation_Stride, [], 100);
F_Umbi_100_Stride = mean(F_Umbi_Interpolation_Stride);
F_Umbi_100_Stride = F_Umbi_100_Stride.';

F_Umbi_Interpolation_Contact = interp1(F_Umbi(1:TO-1), Vector_Interpolation_Contact, 'linear');
F_Umbi_Interpolation_Contact = F_Umbi_Interpolation_Contact';
F_Umbi_Interpolation_Contact = reshape (F_Umbi_Interpolation_Contact, [], 100);
F_Umbi_100_Contact = mean (F_Umbi_Interpolation_Contact);
F_Umbi_100_Contact = F_Umbi_100_Contact.';

F_Umbi_Interpolation_Swing = interp1(F_Umbi(TO:end), Vector_Interpolation_Swing, 'linear');
F_Umbi_Interpolation_Swing = F_Umbi_Interpolation_Swing';
F_Umbi_Interpolation_Swing = reshape (F_Umbi_Interpolation_Swing, [], 100);
F_Umbi_100_Swing = mean (F_Umbi_Interpolation_Swing);
F_Umbi_100_Swing = F_Umbi_100_Swing.';


% DRAG FORCE UMBILICAL NORMALIZED BY BODY  WEIGHT
F_Umbi_WeightNormalized_Interpolation_Stride = interp1(F_Umbi_WeightNormalized, Vector_Interpolation_Stride, 'linear');
F_Umbi_WeightNormalized_Interpolation_Stride = F_Umbi_WeightNormalized_Interpolation_Stride';
F_Umbi_WeightNormalized_Interpolation_Stride = reshape (F_Umbi_WeightNormalized_Interpolation_Stride, [], 100);
F_Umbi_WeightNormalized_100_Stride = mean(F_Umbi_WeightNormalized_Interpolation_Stride);
F_Umbi_WeightNormalized_100_Stride = F_Umbi_WeightNormalized_100_Stride.';

F_Umbi_WeightNormalized_Interpolation_Contact = interp1(F_Umbi_WeightNormalized(1:TO-1), Vector_Interpolation_Contact, 'linear');
F_Umbi_WeightNormalized_Interpolation_Contact = F_Umbi_WeightNormalized_Interpolation_Contact';
F_Umbi_WeightNormalized_Interpolation_Contact = reshape (F_Umbi_WeightNormalized_Interpolation_Contact, [], 100);
F_Umbi_WeightNormalized_100_Contact = mean (F_Umbi_WeightNormalized_Interpolation_Contact);
F_Umbi_WeightNormalized_100_Contact = F_Umbi_WeightNormalized_100_Contact.';

F_Umbi_WeightNormalized_Interpolation_Swing = interp1(F_Umbi_WeightNormalized(TO:end), Vector_Interpolation_Swing, 'linear');
F_Umbi_WeightNormalized_Interpolation_Swing = F_Umbi_WeightNormalized_Interpolation_Swing';
F_Umbi_WeightNormalized_Interpolation_Swing = reshape (F_Umbi_WeightNormalized_Interpolation_Swing, [], 100);
F_Umbi_WeightNormalized_100_Swing = mean (F_Umbi_WeightNormalized_Interpolation_Swing);
F_Umbi_WeightNormalized_100_Swing = F_Umbi_WeightNormalized_100_Swing.';


% DRAG FORCE UMBILICAL NORMALIZED BY SEGMENT FRONTAL AREA
F_Umbi_AreaNormalized_Interpolation_Stride = interp1(F_Umbi_AreaNormalized, Vector_Interpolation_Stride, 'linear');
F_Umbi_AreaNormalized_Interpolation_Stride = F_Umbi_AreaNormalized_Interpolation_Stride';
F_Umbi_AreaNormalized_Interpolation_Stride = reshape (F_Umbi_AreaNormalized_Interpolation_Stride, [], 100);
F_Umbi_AreaNormalized_100_Stride = mean(F_Umbi_AreaNormalized_Interpolation_Stride);
F_Umbi_AreaNormalized_100_Stride = F_Umbi_AreaNormalized_100_Stride.';

F_Umbi_AreaNormalized_Interpolation_Contact = interp1(F_Umbi_AreaNormalized(1:TO-1), Vector_Interpolation_Contact, 'linear');
F_Umbi_AreaNormalized_Interpolation_Contact = F_Umbi_AreaNormalized_Interpolation_Contact';
F_Umbi_AreaNormalized_Interpolation_Contact = reshape (F_Umbi_AreaNormalized_Interpolation_Contact, [], 100);
F_Umbi_AreaNormalized_100_Contact = mean (F_Umbi_AreaNormalized_Interpolation_Contact);
F_Umbi_AreaNormalized_100_Contact = F_Umbi_AreaNormalized_100_Contact.';

F_Umbi_AreaNormalized_Interpolation_Swing = interp1(F_Umbi_AreaNormalized(TO:end), Vector_Interpolation_Swing, 'linear');
F_Umbi_AreaNormalized_Interpolation_Swing = F_Umbi_AreaNormalized_Interpolation_Swing';
F_Umbi_AreaNormalized_Interpolation_Swing = reshape (F_Umbi_AreaNormalized_Interpolation_Swing, [], 100);
F_Umbi_AreaNormalized_100_Swing = mean (F_Umbi_AreaNormalized_Interpolation_Swing);
F_Umbi_AreaNormalized_100_Swing = F_Umbi_AreaNormalized_100_Swing.';


%% SECTION 09
% UPPER LEG (UL) MODEL %

%Ap = ThighC/2(pi)
%Ad = KneeC/2(pi)
%L = Stature - SittingH - TibialH
%DELSH = SittinhH - Stature - TrochantericH
%DELSH: distance from upper segment end to joint center.

%UPPER LEG SIZE MODELING (CYLINDER)
Ap_UL = ThighC/(2*pi);
Ad_UL = KneeC/(2*pi);
L_UL = Stature - SittingH - TibialH;

% DIVISION OF SEGMENT AREA IN dZ PARTS

% dA(z) = 2*((Ad-Ap)z + Ap*L)/L))*dz
% z = (0+dL_UL(i))

%UL length differential
dL_UL=0:0.001:max(L_UL);

%dA_UL = dA(z)

for i=1:length(dL_UL)
    dA_UL(i) =(((Ap_UL-Ad_UL)*(0+dL_UL(i))+(Ad_UL*L_UL))/L_UL)*0.001*2;
end

% TOTAL SEGMENT FRONTAL AREA
A_UL = sum(dA_UL);


% FRONTAL PROJECTED AREA OF EACH dZ PART AREA

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


% LINEAR VELOCITY OF EACH dZ PART

%%v(z)=((z*Vd)+(L-z)*Vp)/L
%UL Velocity
%Vp: Velocity_Throcanter
%Vd: Velocity_Knee
%V_UL = v(z)

% NESTED FOR LOOP: 
% FIRST FOR LOOP IS TO TIME VECTOR
% SECOND FOR LOOP IS TO SEGMENT dL(z)

%V_UL [frames dL_UL]

V_UL = zeros(length(Velocity_Knee), length(dL_UL));   %preallocation: to optimize routine running time

for j=1:1:length(Velocity_Knee)
    for q=1:length(dL_UL)
        V_UL (j,q) = (((0+dL_UL(q))*Velocity_Throcanter(j))+((L_UL-(0+dL_UL(q)))*Velocity_Knee(j)))/L_UL;
    end
end

% DRAG FORCE IN EACH dZ PART

%Cd: drag coefficient;
%W: water density (values from Orselli & Duarte, 2011)
Cd = 1;
W = 1000;

for j=1:size(V_UL,1)
    for q=1:size(V_UL,2)
        dF_UL(j,q) = Cd*W*0.5*(V_UL(j,q).^2)*(dAP_UL(j,q));
    end
end

% TOTAL DRAG FORCE ON SEGMENT: SUM OF DRAG FORCE OF EACH dZ PART
F_UL = (-1)*sum(dF_UL,2);

% FILTERING OF DRAG FORCE
fcut=8;
order=2;
F_UL = matfiltfilt(dt, fcut, order, F_UL);


% TORQUE (MOMENT) DUE DRAG FORCE IN EACH dZ PART
for i=1:size(dF_UL,1)
    for j=1:size(dF_UL,2)
        if Thigh_Angle(i)<=90;
            dT_UL(i,j) = dF_UL(i,j) * (0+dL_UL(j)) * sind(Thigh_Angle(i));
        else
            dT_UL(i,j) = dF_UL(i,j) * (0+dL_UL(j)) * sind(180-Thigh_Angle(i));
       end
   end
end

% TOTAL TORQUE (MOMENT) ON SEGMENT
T_UL = (-1)*sum(dT_UL,2);


% DIVISION OF DRAG FORCE IN CONTACT AND SWING PHASES
F_UL_Contact = F_UL(1:TO-1);
F_UL_Swing = F_UL(TO:length(F_UL));

F_UL_Contact_Max = max(F_UL_Contact);
F_UL_Contact_Min = min(F_UL_Contact);
F_UL_Contact_Mean = mean(F_UL_Contact);

F_UL_Swing_Max = max(F_UL_Swing);
F_UL_Swing_Min = min(F_UL_Swing);
F_UL_Swing_Mean = mean(F_UL_Swing);


% CALCULATION OF RMS DRAG FORCE OF CONTACT AND SWING PHASES
RMS_F_UL_Contact = rms(F_UL_Contact);
RMS_F_UL_Swing = rms(F_UL_Swing);


% NORMALIZATION OF RMS DRAG FORCE BY SEGMENT FRONTAL AREA
RMS_F_UL_AreaNormalized_Contact = RMS_F_UL_Contact/A_UL; 
RMS_F_UL_AreaNormalized_Swing = RMS_F_UL_Swing/A_UL; 

F_UL_AreaNormalized = F_UL/A_UL;

% NORMALIZATION OF DRAG FORCE BY BODY WEIGHT 
F_UL_WeightNormalized = (F_UL/Body_Weight)*100;

F_UL_WeightNormalized_Contact = F_UL_WeightNormalized(1:TO-1);
F_UL_WeightNormalized_Swing = F_UL_WeightNormalized(TO:length(F_UL_WeightNormalized));

RMS_F_UL_WeightNormalized_Contact = rms(F_UL_WeightNormalized_Contact);
RMS_F_UL_WeightNormalized_Swing = rms(F_UL_WeightNormalized_Swing);



    % INTERPOLATION TO 0-100% OF STRIDE DURATION
        
% DRAG FORCE UPPER LEG (UL)
F_UL_Interpolation_Stride = interp1(F_UL, Vector_Interpolation_Stride, 'linear');
F_UL_Interpolation_Stride = F_UL_Interpolation_Stride';
F_UL_Interpolation_Stride = reshape (F_UL_Interpolation_Stride, [], 100);
F_UL_100_Stride = mean(F_UL_Interpolation_Stride);
F_UL_100_Stride = F_UL_100_Stride.';

F_UL_Interpolation_Contact = interp1(F_UL(1:TO-1), Vector_Interpolation_Contact, 'linear');
F_UL_Interpolation_Contact = F_UL_Interpolation_Contact';
F_UL_Interpolation_Contact = reshape (F_UL_Interpolation_Contact, [], 100);
F_UL_100_Contact = mean (F_UL_Interpolation_Contact);
F_UL_100_Contact = F_UL_100_Contact.';

F_UL_Interpolation_Swing = interp1(F_UL(TO:end), Vector_Interpolation_Swing, 'linear');
F_UL_Interpolation_Swing = F_UL_Interpolation_Swing';
F_UL_Interpolation_Swing = reshape (F_UL_Interpolation_Swing, [], 100);
F_UL_100_Swing = mean (F_UL_Interpolation_Swing);
F_UL_100_Swing = F_UL_100_Swing.';


% DRAG FORCE UMBILICAL NORMALIZED BY BODY  WEIGHT
F_UL_WeightNormalized_Interpolation_Stride = interp1(F_UL_WeightNormalized, Vector_Interpolation_Stride, 'linear');
F_UL_WeightNormalized_Interpolation_Stride = F_UL_WeightNormalized_Interpolation_Stride';
F_UL_WeightNormalized_Interpolation_Stride = reshape (F_UL_WeightNormalized_Interpolation_Stride, [], 100);
F_UL_WeightNormalized_100_Stride = mean(F_UL_WeightNormalized_Interpolation_Stride);
F_UL_WeightNormalized_100_Stride = F_UL_WeightNormalized_100_Stride.';

F_UL_WeightNormalized_Interpolation_Contact = interp1(F_UL_WeightNormalized(1:TO-1), Vector_Interpolation_Contact, 'linear');
F_UL_WeightNormalized_Interpolation_Contact = F_UL_WeightNormalized_Interpolation_Contact';
F_UL_WeightNormalized_Interpolation_Contact = reshape (F_UL_WeightNormalized_Interpolation_Contact, [], 100);
F_UL_WeightNormalized_100_Contact = mean (F_UL_WeightNormalized_Interpolation_Contact);
F_UL_WeightNormalized_100_Contact = F_UL_WeightNormalized_100_Contact.';

F_UL_WeightNormalized_Interpolation_Swing = interp1(F_UL_WeightNormalized(TO:end), Vector_Interpolation_Swing, 'linear');
F_UL_WeightNormalized_Interpolation_Swing = F_UL_WeightNormalized_Interpolation_Swing';
F_UL_WeightNormalized_Interpolation_Swing = reshape (F_UL_WeightNormalized_Interpolation_Swing, [], 100);
F_UL_WeightNormalized_100_Swing = mean (F_UL_WeightNormalized_Interpolation_Swing);
F_UL_WeightNormalized_100_Swing = F_UL_WeightNormalized_100_Swing.';


% DRAG FORCE UMBILICAL NORMALIZED BY SEGMENT FRONTAL AREA
F_UL_AreaNormalized_Interpolation_Stride = interp1(F_UL_AreaNormalized, Vector_Interpolation_Stride, 'linear');
F_UL_AreaNormalized_Interpolation_Stride = F_UL_AreaNormalized_Interpolation_Stride';
F_UL_AreaNormalized_Interpolation_Stride = reshape (F_UL_AreaNormalized_Interpolation_Stride, [], 100);
F_UL_AreaNormalized_100_Stride = mean(F_UL_AreaNormalized_Interpolation_Stride);
F_UL_AreaNormalized_100_Stride = F_UL_AreaNormalized_100_Stride.';

F_UL_AreaNormalized_Interpolation_Contact = interp1(F_UL_AreaNormalized(1:TO-1), Vector_Interpolation_Contact, 'linear');
F_UL_AreaNormalized_Interpolation_Contact = F_UL_AreaNormalized_Interpolation_Contact';
F_UL_AreaNormalized_Interpolation_Contact = reshape (F_UL_AreaNormalized_Interpolation_Contact, [], 100);
F_UL_AreaNormalized_100_Contact = mean (F_UL_AreaNormalized_Interpolation_Contact);
F_UL_AreaNormalized_100_Contact = F_UL_AreaNormalized_100_Contact.';

F_UL_AreaNormalized_Interpolation_Swing = interp1(F_UL_AreaNormalized(TO:end), Vector_Interpolation_Swing, 'linear');
F_UL_AreaNormalized_Interpolation_Swing = F_UL_AreaNormalized_Interpolation_Swing';
F_UL_AreaNormalized_Interpolation_Swing = reshape (F_UL_AreaNormalized_Interpolation_Swing, [], 100);
F_UL_AreaNormalized_100_Swing = mean (F_UL_AreaNormalized_Interpolation_Swing);
F_UL_AreaNormalized_100_Swing = F_UL_AreaNormalized_100_Swing.';

%% SECTION 10
% LOWER LEG (LL) MODEL %

%Ap = KneeC/2(pi)
%Ad = AnkleC/2(pi)
%SL = TibialH - SphrH

% LOWER LEG SIZE MODELING (CYLINDER)
Ap_LL = KneeC/(2*pi);
Ad_LL = AnkleC/(2*pi);
L_LL = TibialH - SphrH;

% DIVISION OF SEGMENT AREA IN dZ PARTS

% dA(z) = 2*((Ad-Ap)z + Ap*L)/L))*dz
% z = (0+dL_LL(i))

%LL length differential
dL_LL=0:0.001:max(L_LL);

%dA_LL = dA(z)
for i=1:length(dL_LL)
    dA_LL(i) =(((Ap_LL-Ad_LL)*(0+dL_LL(i))+(Ad_LL*L_LL))/L_LL)*0.001*2;
end


% TOTAL SEGMENT FRONTAL AREA
A_LL = sum(dA_LL);


% FRONTAL PROJECTED AREA OF EACH dZ PART AREA

% dAP_LL [frames dA_LLparts] = dAP(z)
for i=1:length(Shank_Angle);
    
    for j=1:length(dA_LL);
        if Shank_Angle(i)<=90;
            dAP_LL(i,j)=dA_LL(j)*sind(Shank_Angle(i));
        else
            dAP_LL(i,j)=dA_LL(j)*sind(180-Shank_Angle(i));
        end
    end
end


% LINEAR VELOCITY OF EACH dZ PART

% v(z)=((z*Vd)+(L-z)*Vp)/L
% LL Velocity
% Vp: Velocity_Knee
% Vd: Velocity_Ankle
% V_LL = v(z)

% NESTED FOR LOOP: 
% FIRST FOR LOOP IS TO TIME VECTOR
% SECOND FOR LOOP IS TO SEGMENT dL(z)

%V_LL [frames dL_LL]

V_LL = zeros(length(Velocity_Knee), length(dL_LL));   %preallocation: to optimize routine running time

for j=1:1:length(Velocity_Knee)
    
    for q=1:length(dL_LL)
        V_LL (j,q) = (((0+dL_LL(q))*Velocity_Ankle(j))+((L_LL-(0+dL_LL(q)))*Velocity_Knee(j)))/L_LL;
    end
end

% DRAG FORCE IN EACH dZ PART

%Cd: drag coefficient;
%W: water density (values from Orselli & Duarte, 2011)
Cd = 1;
W = 1000;

for j=1:size(V_LL,1)
    for q=1:size(V_LL,2)
        dF_LL(j,q) = Cd*W*0.5*(V_LL(j,q).^2)*(dAP_LL(j,q));
    end
end


% TOTAL DRAG FORCE ON SEGMENT: SUM OF DRAG FORCE OF EACH dZ PART
F_LL = (-1)*sum(dF_LL,2);

% FILTERING OF DRAG FORCE
fcut=8;
order=2;
F_LL = matfiltfilt(dt, fcut, order, F_LL);


% TORQUE (MOMENT) DUE DRAG FORCE IN EACH dZ PART
for i=1:size(dF_LL,1)
    for j=1:size(dF_LL,2)
        if Shank_Angle(i)<=90;
            dT_LL(i,j) = dF_LL(i,j) * (0+dL_LL(j)) * sind(Shank_Angle(i));
        else
            dT_LL(i,j) = dF_LL(i,j) * (0+dL_LL(j)) * sind(180-Shank_Angle(i));
       end
   end
end

% TOTAL TORQUE (MOMENT) ON SEGMENT
T_LL = (-1)*sum(dT_LL,2);


% DIVISION OF DRAG FORCE IN CONTACT AND SWING PHASES
F_LL_Contact = F_LL(1:TO-1);
F_LL_Swing = F_LL(TO:length(F_LL));

F_LL_Contact_Max = max(F_LL_Contact);
F_LL_Contact_Min = min(F_LL_Contact);
F_LL_Contact_Mean = mean(F_LL_Contact);

F_LL_Swing_Max = max(F_LL_Swing);
F_LL_Swing_Min = min(F_LL_Swing);
F_LL_Swing_Mean = mean(F_LL_Swing);


% CALCULATION OF RMS DRAG FORCE OF CONTACT AND SWING PHASES
RMS_F_LL_Contact = rms(F_LL_Contact);
RMS_F_LL_Swing = rms(F_LL_Swing);


% NORMALIZATION OF RMS DRAG FORCE BY SEGMENT FRONTAL AREA
RMS_F_LL_AreaNormalized_Contact = RMS_F_LL_Contact/A_LL; 
RMS_F_LL_AreaNormalized_Swing = RMS_F_LL_Swing/A_LL; 

F_LL_AreaNormalized = F_LL/A_LL;

% NORMALIZATION OF DRAG FORCE BY BODY WEIGHT 
F_LL_WeightNormalized = (F_LL/Body_Weight)*100;

F_LL_WeightNormalized_Contact = F_LL_WeightNormalized(1:TO-1);
F_LL_WeightNormalized_Swing = F_LL_WeightNormalized(TO:length(F_LL_WeightNormalized));

RMS_F_LL_WeightNormalized_Contact = rms(F_LL_WeightNormalized_Contact);
RMS_F_LL_WeightNormalized_Swing = rms(F_LL_WeightNormalized_Swing);


     % INTERPOLATION TO 0-100% OF STRIDE DURATION
        
% DRAG FORCE LOWER LEG (LL)
F_LL_Interpolation_Stride = interp1(F_LL, Vector_Interpolation_Stride, 'linear');
F_LL_Interpolation_Stride = F_LL_Interpolation_Stride';
F_LL_Interpolation_Stride = reshape (F_LL_Interpolation_Stride, [], 100);
F_LL_100_Stride = mean(F_LL_Interpolation_Stride);
F_LL_100_Stride = F_LL_100_Stride.';

F_LL_Interpolation_Contact = interp1(F_LL(1:TO-1), Vector_Interpolation_Contact, 'linear');
F_LL_Interpolation_Contact = F_LL_Interpolation_Contact';
F_LL_Interpolation_Contact = reshape (F_LL_Interpolation_Contact, [], 100);
F_LL_100_Contact = mean (F_LL_Interpolation_Contact);
F_LL_100_Contact = F_LL_100_Contact.';

F_LL_Interpolation_Swing = interp1(F_LL(TO:end), Vector_Interpolation_Swing, 'linear');
F_LL_Interpolation_Swing = F_LL_Interpolation_Swing';
F_LL_Interpolation_Swing = reshape (F_LL_Interpolation_Swing, [], 100);
F_LL_100_Swing = mean (F_LL_Interpolation_Swing);
F_LL_100_Swing = F_LL_100_Swing.';


% DRAG FORCE LOWER LEG NORMALIZED BY BODY  WEIGHT
F_LL_WeightNormalized_Interpolation_Stride = interp1(F_LL_WeightNormalized, Vector_Interpolation_Stride, 'linear');
F_LL_WeightNormalized_Interpolation_Stride = F_LL_WeightNormalized_Interpolation_Stride';
F_LL_WeightNormalized_Interpolation_Stride = reshape (F_LL_WeightNormalized_Interpolation_Stride, [], 100);
F_LL_WeightNormalized_100_Stride = mean(F_LL_WeightNormalized_Interpolation_Stride);
F_LL_WeightNormalized_100_Stride = F_LL_WeightNormalized_100_Stride.';

F_LL_WeightNormalized_Interpolation_Contact = interp1(F_LL_WeightNormalized(1:TO-1), Vector_Interpolation_Contact, 'linear');
F_LL_WeightNormalized_Interpolation_Contact = F_LL_WeightNormalized_Interpolation_Contact';
F_LL_WeightNormalized_Interpolation_Contact = reshape (F_LL_WeightNormalized_Interpolation_Contact, [], 100);
F_LL_WeightNormalized_100_Contact = mean (F_LL_WeightNormalized_Interpolation_Contact);
F_LL_WeightNormalized_100_Contact = F_LL_WeightNormalized_100_Contact.';

F_LL_WeightNormalized_Interpolation_Swing = interp1(F_LL_WeightNormalized(TO:end), Vector_Interpolation_Swing, 'linear');
F_LL_WeightNormalized_Interpolation_Swing = F_LL_WeightNormalized_Interpolation_Swing';
F_LL_WeightNormalized_Interpolation_Swing = reshape (F_LL_WeightNormalized_Interpolation_Swing, [], 100);
F_LL_WeightNormalized_100_Swing = mean (F_LL_WeightNormalized_Interpolation_Swing);
F_LL_WeightNormalized_100_Swing = F_LL_WeightNormalized_100_Swing.';


% DRAG FORCE LOWER LEG NORMALIZED BY SEGMENT FRONTAL AREA
F_LL_AreaNormalized_Interpolation_Stride = interp1(F_LL_AreaNormalized, Vector_Interpolation_Stride, 'linear');
F_LL_AreaNormalized_Interpolation_Stride = F_LL_AreaNormalized_Interpolation_Stride';
F_LL_AreaNormalized_Interpolation_Stride = reshape (F_LL_AreaNormalized_Interpolation_Stride, [], 100);
F_LL_AreaNormalized_100_Stride = mean(F_LL_AreaNormalized_Interpolation_Stride);
F_LL_AreaNormalized_100_Stride = F_LL_AreaNormalized_100_Stride.';

F_LL_AreaNormalized_Interpolation_Contact = interp1(F_LL_AreaNormalized(1:TO-1), Vector_Interpolation_Contact, 'linear');
F_LL_AreaNormalized_Interpolation_Contact = F_LL_AreaNormalized_Interpolation_Contact';
F_LL_AreaNormalized_Interpolation_Contact = reshape (F_LL_AreaNormalized_Interpolation_Contact, [], 100);
F_LL_AreaNormalized_100_Contact = mean (F_LL_AreaNormalized_Interpolation_Contact);
F_LL_AreaNormalized_100_Contact = F_LL_AreaNormalized_100_Contact.';

F_LL_AreaNormalized_Interpolation_Swing = interp1(F_LL_AreaNormalized(TO:end), Vector_Interpolation_Swing, 'linear');
F_LL_AreaNormalized_Interpolation_Swing = F_LL_AreaNormalized_Interpolation_Swing';
F_LL_AreaNormalized_Interpolation_Swing = reshape (F_LL_AreaNormalized_Interpolation_Swing, [], 100);
F_LL_AreaNormalized_100_Swing = mean (F_LL_AreaNormalized_Interpolation_Swing);
F_LL_AreaNormalized_100_Swing = F_LL_AreaNormalized_100_Swing.';

%% SECTION 11
% FOOT (FT) MODEL %

%Ap = SphrH*0.5
%Ad = ((1.287*L)-L)/((1.287*L)/Ap_Ft)           %Calculate from trigonometry (Triangles Similarity), CG from Foot is at Ap/3 from axis, and at 0.429L from Proximal Joint.
%SL = FootL

% FOOT SIZE MODELING (CYLINDER)
Ap_FT = SphrH*0.5;
L_FT = FootL;
Ad_FT = ((1.287*L_FT)-L_FT)/((1.287*L_FT)/Ap_FT);

% DIVISION OF SEGMENT AREA IN dZ PARTS

%%dA(z) = 2*((Ad-Ap)z + Ap*L)/L))*dz
% z = (0+dL_FT(i))

%FT length differential
dL_FT=0:0.001:max(L_FT);

%dA_FT = dA(z)
for i=1:length(dL_FT)
    dA_FT(i) =(((Ap_FT-Ad_FT)*(0+dL_FT(i))+(Ad_FT*L_FT))/L_FT)*0.001*2;
end


% TOTAL SEGMENT FRONTAL AREA
A_FT = sum(dA_FT);


% FRONTAL PROJECTED AREA OF EACH dZ PART AREA

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


% LINEAR VELOCITY OF EACH dZ PART

% v(z)=((z*Vd)+(L-z)*Vp)/L
% FT Velocity
% Vp: Velocity_Knee
% Vd: Velocity_Toe
% V_FT = v(z)

% NESTED FOR LOOP: 
% FIRST FOR LOOP IS TO TIME VECTOR
% SECOND FOR LOOP IS TO SEGMENT dL(z)

%V_FT [frames dL_FT]

V_FT = zeros(length(Velocity_Ankle), length(dL_FT));   %preallocation: to optimize routine running time

for j=1:1:length(Velocity_Ankle)
    
    for q=1:length(dL_FT)
        V_FT (j,q) = (((0+dL_FT(q))*Velocity_Toe(j))+((L_FT-(0+dL_FT(q)))*Velocity_Ankle(j)))/L_FT;
    end
end


% DRAG FORCE IN EACH dZ PART

%Cd: drag coefficient;
%W: water density (values from Orselli & Duarte, 2011)
Cd = 1;
W = 1000;

for j=1:size(V_FT,1)
    for q=1:size(V_FT,2)
        dF_FT(j,q) = Cd*W*0.5*(V_FT(j,q).^2)*(dAP_FT(j,q));
    end
end


% TOTAL DRAG FORCE ON SEGMENT: SUM OF DRAG FORCE OF EACH dZ PART
F_FT = (-1)*sum(dF_FT,2);


% FILTERING OF DRAG FORCE
fcut=8;
order=2;
F_FT = matfiltfilt(dt, fcut, order, F_FT);

% TORQUE (MOMENT) DUE DRAG FORCE IN EACH dZ PART
for i=1:size(dF_FT,1)
    for j=1:size(dF_FT,2)
        if Foot_Angle(i)<=90;
            dT_FT(i,j) = dF_FT(i,j) * (0+dL_FT(j)) * sind(Foot_Angle(i));
        else
            dT_FT(i,j) = dF_FT(i,j) * (0+dL_FT(j)) * sind(180-Foot_Angle(i));
       end
   end
end

% TOTAL TORQUE (MOMENT) ON SEGMENT
T_FT = (-1)*sum(dT_FT,2);



% DIVISION OF DRAG FORCE IN CONTACT AND SWING PHASES
F_FT_Contact = F_FT(1:TO-1);
F_FT_100_Swing = F_FT(TO:length(F_FT));

F_FT_Contact_Max = max(F_FT_Contact);
F_FT_Contact_Min = min(F_FT_Contact);
F_FT_Contact_Mean = mean(F_FT_Contact);

F_FT_Swing_Max = max(F_FT_100_Swing);
F_FT_Swing_Min = min(F_FT_100_Swing);
F_FT_Swing_Mean = mean(F_FT_100_Swing);


% CALCULATION OF RMS DRAG FORCE OF CONTACT AND SWING PHASES
RMS_F_FT_Contact = rms(F_FT_Contact);
RMS_F_FT_Swing = rms(F_FT_100_Swing);


% NORMALIZATION OF RMS DRAG FORCE BY SEGMENT FRONTAL AREA
RMS_F_FT_AreaNormalized_Contact = RMS_F_FT_Contact/A_FT; 
RMS_F_FT_AreaNormalized_Swing = RMS_F_FT_Swing/A_FT; 

F_FT_AreaNormalized = F_FT/A_FT;

% NORMALIZATION OF DRAG FORCE BY BODY WEIGHT 
F_FT_WeightNormalized = (F_FT/Body_Weight)*100;

F_FT_WeightNormalized_Contact = F_FT_WeightNormalized(1:TO-1);
F_FT_WeightNormalized_Swing = F_FT_WeightNormalized(TO:length(F_FT_WeightNormalized));

RMS_F_FT_WeightNormalized_Contact = rms(F_FT_WeightNormalized_Contact);
RMS_F_FT_WeightNormalized_Swing = rms(F_FT_WeightNormalized_Swing);


    % INTERPOLATION TO 0-100% OF STRIDE DURATION
        
% DRAG FORCE FOOT (FT)
F_FT_Interpolation_Stride = interp1(F_FT, Vector_Interpolation_Stride, 'linear');
F_FT_Interpolation_Stride = F_FT_Interpolation_Stride';
F_FT_Interpolation_Stride = reshape (F_FT_Interpolation_Stride, [], 100);
F_FT_100_Stride = mean(F_FT_Interpolation_Stride);
F_FT_100_Stride = F_FT_100_Stride.';

F_FT_Interpolation_Contact = interp1(F_FT(1:TO-1), Vector_Interpolation_Contact, 'linear');
F_FT_Interpolation_Contact = F_FT_Interpolation_Contact';
F_FT_Interpolation_Contact = reshape (F_FT_Interpolation_Contact, [], 100);
F_FT_100_Contact = mean (F_FT_Interpolation_Contact);
F_FT_100_Contact = F_FT_100_Contact.';

F_FT_Interpolation_Swing = interp1(F_FT(TO:end), Vector_Interpolation_Swing, 'linear');
F_FT_Interpolation_Swing = F_FT_Interpolation_Swing';
F_FT_Interpolation_Swing = reshape (F_FT_Interpolation_Swing, [], 100);
F_FT_100_Swing = mean (F_FT_Interpolation_Swing);
F_FT_100_Swing = F_FT_100_Swing.';


% DRAG FORCE FOOT NORMALIZED BY BODY  WEIGHT
F_FT_WeightNormalized_Interpolation_Stride = interp1(F_FT_WeightNormalized, Vector_Interpolation_Stride, 'linear');
F_FT_WeightNormalized_Interpolation_Stride = F_FT_WeightNormalized_Interpolation_Stride';
F_FT_WeightNormalized_Interpolation_Stride = reshape (F_FT_WeightNormalized_Interpolation_Stride, [], 100);
F_FT_WeightNormalized_100_Stride = mean(F_FT_WeightNormalized_Interpolation_Stride);
F_FT_WeightNormalized_100_Stride = F_FT_WeightNormalized_100_Stride.';

F_FT_WeightNormalized_Interpolation_Contact = interp1(F_FT_WeightNormalized(1:TO-1), Vector_Interpolation_Contact, 'linear');
F_FT_WeightNormalized_Interpolation_Contact = F_FT_WeightNormalized_Interpolation_Contact';
F_FT_WeightNormalized_Interpolation_Contact = reshape (F_FT_WeightNormalized_Interpolation_Contact, [], 100);
F_FT_WeightNormalized_100_Contact = mean (F_FT_WeightNormalized_Interpolation_Contact);
F_FT_WeightNormalized_100_Contact = F_FT_WeightNormalized_100_Contact.';

F_FT_WeightNormalized_Interpolation_Swing = interp1(F_FT_WeightNormalized(TO:end), Vector_Interpolation_Swing, 'linear');
F_FT_WeightNormalized_Interpolation_Swing = F_FT_WeightNormalized_Interpolation_Swing';
F_FT_WeightNormalized_Interpolation_Swing = reshape (F_FT_WeightNormalized_Interpolation_Swing, [], 100);
F_FT_WeightNormalized_100_Swing = mean (F_FT_WeightNormalized_Interpolation_Swing);
F_FT_WeightNormalized_100_Swing = F_FT_WeightNormalized_100_Swing.';


% DRAG FORCE FOOT NORMALIZED BY SEGMENT FRONTAL AREA
F_FT_AreaNormalized_Interpolation_Stride = interp1(F_FT_AreaNormalized, Vector_Interpolation_Stride, 'linear');
F_FT_AreaNormalized_Interpolation_Stride = F_FT_AreaNormalized_Interpolation_Stride';
F_FT_AreaNormalized_Interpolation_Stride = reshape (F_FT_AreaNormalized_Interpolation_Stride, [], 100);
F_FT_AreaNormalized_100_Stride = mean(F_FT_AreaNormalized_Interpolation_Stride);
F_FT_AreaNormalized_100_Stride = F_FT_AreaNormalized_100_Stride.';

F_FT_AreaNormalized_Interpolation_Contact = interp1(F_FT_AreaNormalized(1:TO-1), Vector_Interpolation_Contact, 'linear');
F_FT_AreaNormalized_Interpolation_Contact = F_FT_AreaNormalized_Interpolation_Contact';
F_FT_AreaNormalized_Interpolation_Contact = reshape (F_FT_AreaNormalized_Interpolation_Contact, [], 100);
F_FT_AreaNormalized_100_Contact = mean (F_FT_AreaNormalized_Interpolation_Contact);
F_FT_AreaNormalized_100_Contact = F_FT_AreaNormalized_100_Contact.';

F_FT_AreaNormalized_Interpolation_Swing = interp1(F_FT_AreaNormalized(TO:end), Vector_Interpolation_Swing, 'linear');
F_FT_AreaNormalized_Interpolation_Swing = F_FT_AreaNormalized_Interpolation_Swing';
F_FT_AreaNormalized_Interpolation_Swing = reshape (F_FT_AreaNormalized_Interpolation_Swing, [], 100);
F_FT_AreaNormalized_100_Swing = mean (F_FT_AreaNormalized_Interpolation_Swing);
F_FT_AreaNormalized_100_Swing = F_FT_AreaNormalized_100_Swing.';

%% SECTION 12
% DRAG FORCES SUMMATION

F_Contact_Total = RMS_F_Xiph_Contact + RMS_F_UL_Contact + RMS_F_LL_Contact + RMS_F_FT_Contact;
F_Swing_Total = RMS_F_Xiph_Swing + RMS_F_UL_Swing + RMS_F_LL_Swing + RMS_F_FT_Swing;
F_Stride_Total = F_Contact_Total + F_Swing_Total; 

F_Contact_Total_WeightNormalized = (F_Contact_Total/Body_Weight)*100;
F_Swing_Total_WeightNormalized = (F_Swing_Total/Body_Weight)*100;
F_Stride_Total_WeightNormalized = (F_Stride_Total/Body_Weight)*100;
%%  SECTION 13
% GROUND REACTION FORCES ESTIMATION %
% (Haupenthal et al., 2019)

GRFV_PEAK = 566.682 - 922.059*(Stature/Immersion_depth) + 3.321 * (Mass) + (176.369 * Speed_Mean);
GRFAP_PEAK = -164.160 + (243.626* Speed_Mean) + (2.443*((ThighC + KneeC)/2));

GRFV_PEAK_WeightNormalized = (GRFV_PEAK/Body_Weight)*100;
GRFAP_PEAK_WeightNormalized = (GRFAP_PEAK/Body_Weight)*100;

%% SECTION 14
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

%(Hreljac & Marshal, 2000)

%% SECTION 15 
% EXPORT DATA %

               %%% SPATIOTEMPORAL DATA        
Export_File_Path_Full_Spatiotemporal = [Export_File_Path '\' File_Name '_SPATIOTEMPORAL' '.xls'];
      
Export_Data_SpatioTemporal = [Speed_Mean Stride_Frequency Stride_Lenght Stride_Duration Contact_Phase_Duration Duty_Factor TO_in_Percentage];
Export_Header_Spatiotemporal = {'Speed_Mean (m/s)', 'Stride_Frequency (Hz)', 'Stride_Lenght (m)', 'Stride_Duration (s)', 'Contact_Phase_Duration (s)', 'Duty Factor (%)', 'TO in % of Stride'};

Export_Excel_Spatiotemporal = [Export_Header_Spatiotemporal; num2cell(Export_Data_SpatioTemporal)];
xlswrite (Export_File_Path_Full_Spatiotemporal, Export_Excel_Spatiotemporal); 
    

            %%% ANGULAR DATA
Export_File_Path_Full_Angular = [Export_File_Path '\' File_Name '_ANGULAR' '.xls'];
      
Export_Data_Angular = [Foot_Angle_ROM Shank_Angle_ROM Thigh_Angle_ROM Trunk_Angle_ROM Foot_Angular_Speed_Max Shank_Angular_Speed_Max Thigh_Angular_Speed_Max Trunk_Angular_Speed_Max Ankle_Joint_Angle_ROM Knee_Joint_Angle_ROM Hip_Joint_Angle_ROM Ankle_Joint_Angular_Speed_Max Knee_Joint_Angular_Speed_Max Hip_Joint_Angular_Speed_Max];
Export_Header_Angular = {'Foot_Angle_ROM (º)', 'Shank_Angle_ROM (º)', 'Thigh_Angle_ROM (º)', 'Trunk_Angle_ROM (º)', 'Foot_Angular_Speed_Max (º/s)', 'Shank_Angular_Speed_Max (º/s)', 'Thigh_Angular_Speed_Max (º/s)', 'Trunk_Angular_Speed_Max (º/s)', 'Ankle_Joint_Angle_ROM (º)', 'Knee_Joint_Angle_ROM (º)', 'Hip_Joint_Angle_ROM (º)', 'Ankle_Joint_Angular_Speed_Max (º/s)', 'Knee_Joint_Angular_Speed_Max (º/s)', 'Hip_Joint_Angular_Speed_Max (º/s)'};

Export_Excel_Angular = [Export_Header_Angular; num2cell(Export_Data_Angular)];
xlswrite (Export_File_Path_Full_Angular, Export_Excel_Angular); 

% COORDINATION ANGULAR

Export_File_Path_Full_Coordination = [Export_File_Path '\' File_Name '_COORDINATION_ANGULAR' '.xls'];
      
Export_Data_Coordination = [Relative_Phase_Foot_Shank Relative_Phase_Shank_Thigh Coupling_Angle_Shank_Foot_Stride Coupling_Angle_Shank_Foot_Contact Coupling_Angle_Shank_Foot_Swing Coupling_Angle_Thigh_Shank_Stride Coupling_Angle_Thigh_Shank_Contact Coupling_Angle_Thigh_Shank_Swing];
Export_Header_Coordination = {'Relative_Phase_Foot_Shank', 'Relative_Phase_Shank_Thigh', 'Coupling_Angle_Shank_Foot_Stride', 'Coupling_Angle_Shank_Foot_Contact', 'Coupling_Angle_Shank_Foot_Swing', 'Coupling_Angle_Thigh_Shank_Stride', 'Coupling_Angle_Thigh_Shank_Contact', 'Coupling_Angle_Thigh_Shank_Swing'};

Export_Excel_Coordination = [Export_Header_Coordination; num2cell(Export_Data_Coordination)];
xlswrite (Export_File_Path_Full_Coordination, Export_Excel_Coordination); 


            %%% KINETIC DATA: DRAG FORCES AND GRF ESTIMATION
    
% TOTAL KINETIC
Output_File_Path_Full_Kinetic_Total = [Export_File_Path '\' File_Name '_KINETIC_TOTAL' '.xls'];
      
Export_Data_Kinetic_Total= [F_Stride_Total F_Contact_Total F_Swing_Total F_Stride_Total_WeightNormalized F_Contact_Total_WeightNormalized F_Swing_Total_WeightNormalized GRFV_PEAK GRFAP_PEAK GRFV_PEAK_WeightNormalized GRFAP_PEAK_WeightNormalized];
Export_Header_Kinetic_Total = {'F_Stride_Total (N)', 'F_Contact_Total (N)', 'F_Swing_Total (N)','F_Stride_Total_Weight_Normalized (%N/BW)', 'F_Contact_Total_Weight_Normalized (%N/BW)', 'F_Swing_Total_Weight_Normalized (%N/BW)' ,'GRFV_PEAK (N)', 'GRAP_PEAK (N)', 'GRFV_PEAK_NORMALIZED(%BW/N)', 'GRFAP_PEAK_NORMALIZED(%BW/N)'};

Export_Excel_Kinetic_Total = [Export_Header_Kinetic_Total; num2cell(Export_Data_Kinetic_Total)];
xlswrite (Output_File_Path_Full_Kinetic_Total, Export_Excel_Kinetic_Total); 


% DETAILED KINETIC
Output_File_Path_Full_Kinetic_Detailed = [Export_File_Path '\' File_Name '_KINETIC_DETAILED' '.xls'];
      
Export_Data_Kinetic_Detailed= [RMS_F_Xiph_Contact RMS_F_Xiph_Swing RMS_F_Umbi_Contact RMS_F_Umbi_Swing RMS_F_UL_Contact RMS_F_UL_Swing RMS_F_LL_Contact RMS_F_LL_Swing RMS_F_FT_Contact RMS_F_FT_Swing RMS_F_Xiph_AreaNormalized_Contact RMS_F_Xiph_AreaNormalized_Swing RMS_F_Umbi_AreaNormalized_Contact RMS_F_Umbi_AreaNormalized_Swing RMS_F_UL_AreaNormalized_Contact RMS_F_UL_AreaNormalized_Swing RMS_F_LL_AreaNormalized_Contact RMS_F_LL_AreaNormalized_Swing RMS_F_FT_AreaNormalized_Contact RMS_F_FT_AreaNormalized_Swing RMS_F_Xiph_WeightNormalized_Contact RMS_F_Xiph_WeightNormalized_Swing RMS_F_Umbi_WeightNormalized_Contact RMS_F_Umbi_WeightNormalized_Swing RMS_F_UL_WeightNormalized_Contact RMS_F_UL_WeightNormalized_Swing RMS_F_LL_WeightNormalized_Contact RMS_F_LL_WeightNormalized_Swing RMS_F_FT_WeightNormalized_Contact RMS_F_FT_WeightNormalized_Swing];
Export_Header_Kinetic_Detailed = {'RMS_F_Xiph_Contact (N)','RMS_F_Xiph_Swing (N)', 'RMS_F_Umbi_Contact (N)', 'RMS_F_Umbi_Swing (N)' , 'RMS_F_UL_Contact (N)','RMS_F_UL_Swing (N)','RMS_F_LL_Contact (N)', 'RMS_F_LL_Swing (N)', 'RMS_F_FT_Contact (N)', 'RMS_F_FT_Swing (N)','RMS_F_Xiph_AreaNormalized_Contact(N/m2)', 'RMS_F_Xiph_AreaNormalized_Swing(N/m2)', 'RMS_F_Umbi_AreaNormalized_Contact(N/m2)', 'RMS_F_Umbi_AreaNormalized_Swing(N/m2)', 'RMS_F_UL_AreaNormalized_Contact(N/m2)', 'RMS_F_UL_AreaNormalized_Swing(N/m2)', 'RMS_F_LL_AreaNormalized_Contact(N/m2)', 'RMS_F_LL_AreaNormalized_Swing(N/m2)', 'RMS_F_FT_AreaNormalized_Contact(N/m2)', 'RMS_F_FT_AreaNormalized_Swing(N/m2)', 'RMS_F_Xiph_WeightNormalized_Contact(%N/BW)', 'RMS_F_Xiph_WeightNormalized_Swing(%N/BW)', 'RMS_F_Umbi_WeightNormalized_Contact(%N/BW)', 'RMS_F_Umbi_WeightNormalized_Swing(%N/BW)', 'RMS_F_UL_WeightNormalized_Contact(%N/BW)', 'RMS_F_UL_WeightNormalized_Swing(%N/BW)', 'RMS_F_LL_WeightNormalized_Contact(%N/BW)', 'RMS_F_LL_WeightNormalized_Swing(%N/BW)', 'RMS_F_FT_WeightNormalized_Contact(%N/BW)', 'RMS_F_FT_WeightNormalized_Swing(%N/BW)'};

Export_Excel_Kinetic_Detailed = [Export_Header_Kinetic_Detailed; num2cell(Export_Data_Kinetic_Detailed)];
xlswrite (Output_File_Path_Full_Kinetic_Detailed, Export_Excel_Kinetic_Detailed);


            %%% CURVES
            
    % CURVES ANGULAR 0-100%    
Output_File_Path_Full_Curves100_Angular = [Export_File_Path '\' File_Name '_CURVES100_ANGULAR' '.xls'];
      
Data_Curves100_Angular = [Foot_Angle_100_Stride Foot_Angle_100_Contact Foot_Angle_100_Swing Shank_Angle_100_Stride Shank_Angle_100_Contact Shank_Angle_100_Swing Thigh_Angle_100_Stride Thigh_Angle_100_Contact Thigh_Angle_100_Swing Ankle_Joint_Angle_100_Stride Ankle_Joint_Angle_100_Contact Ankle_Joint_Angle_100_Swing Knee_Joint_Angle_100_Stride Knee_Joint_Angle_100_Contact Knee_Joint_Angle_100_Swing Hip_Joint_Angle_100_Stride Hip_Joint_Angle_100_Contact Hip_Joint_Angle_100_Swing];
Header_Curves100_Angular = {'Foot_Angle_100_Stride (º)', 'Foot_Angle_100_Contact (º)', 'Foot_Angle_100_Swing (º)', 'Shank_Angle_100_Stride (º)', 'Shank_Angle_100_Contact (º)', 'Shank_Angle_100_Swing (º)', 'Thigh_Angle_100_Stride (º)', 'Thigh_Angle_100_Contact (º)', 'Thigh_Angle_100_Swing (º)', 'Ankle_Joint_Angle_100_Stride (º)', 'Ankle_Joint_Angle_100_Contact (º)', 'Ankle_Joint_Angle_100_Swing (º)', 'Knee_Joint_Angle_100_Stride (º)', 'Knee_Joint_Angle_100_Contact (º)', 'Knee_Joint_Angle_100_Swing (º)', 'Hip_Joint_Angle_100_Stride (º)', 'Hip_Joint_Angle_100_Contact (º)','Hip_Joint_Angle_100_Swing (º)'};

Export_Excel_Curves100_Angular = [Header_Curves100_Angular; num2cell(Data_Curves100_Angular)];
xlswrite (Output_File_Path_Full_Curves100_Angular, Export_Excel_Curves100_Angular); 



    % CURVES KINETIC 0-100% 
    
% ABSOLUTE CURVES    
Output_File_Path_Full_Curves100_Kinetic_Absolute = [Export_File_Path '\' File_Name '_CURVES100_KINETIC_ABSOLUTE' '.xls'];
      
Data_Curves100_Kinetic_Absolute = [F_Xiph_100_Stride F_Xiph_100_Contact F_Xiph_100_Swing F_Umbi_100_Stride F_Umbi_100_Contact F_Umbi_100_Swing F_UL_100_Stride F_UL_100_Contact F_UL_100_Swing F_LL_100_Stride F_LL_100_Contact F_LL_100_Swing F_FT_100_Stride F_FT_100_Contact F_FT_100_Swing ];
Header_Curves100_Kinetic_Absolute = {'F_Xiph_100_Stride (N)', 'F_Xiph_100_Contact (N)', 'F_Xiph_100_Swing (N)', 'F_Umbi_100_Stride (N)', 'F_Umbi_100_Contact (N)', 'F_Umbi_100_Swing (N)', 'F_UL_100_Stride (N)', 'F_UL_100_Contact (N)', 'F_UL_100_Swing (N)', 'F_LL_100_Stride (N)', 'F_LL_100_Contact (N)', 'F_LL_100_Swing (N)', 'F_FT_100_Stride (N)', 'F_FT_100_Contact (N)', 'F_FT_100_Swing (N)'};

Export_Excel_Curves100_Kinetic_Absolute = [Header_Curves100_Kinetic_Absolute; num2cell(Data_Curves100_Kinetic_Absolute)];
xlswrite (Output_File_Path_Full_Curves100_Kinetic_Absolute, Export_Excel_Curves100_Kinetic_Absolute); 


% NORMALIZED CURVES BY BODY WEIGHT  
Output_File_Path_Full_Curves100_Kinetic_Normalized_Weight = [Export_File_Path '\' File_Name '_CURVES100_KINETIC_WEIGHT_NORMALIZED' '.xls'];
      
Data_Curves100_Kinetic_Normalized_Weight = [F_Xiph_WeightNormalized_100_Stride F_Xiph_WeightNormalized_100_Contact F_Xiph_WeightNormalized_100_Swing F_Umbi_WeightNormalized_100_Stride F_Umbi_WeightNormalized_100_Contact F_Umbi_WeightNormalized_100_Swing F_UL_WeightNormalized_100_Stride F_UL_WeightNormalized_100_Contact F_UL_WeightNormalized_100_Swing F_LL_WeightNormalized_100_Stride F_LL_WeightNormalized_100_Contact F_LL_WeightNormalized_100_Swing F_FT_WeightNormalized_100_Stride F_FT_WeightNormalized_100_Contact F_FT_WeightNormalized_100_Swing];
Header_Curves100_Kinetic_Normalized_Weight = {'F_Xiph_WeightNormalized_100_Stride(%N/BW)', 'F_Xiph_WeightNormalized_100_Contact(%N/BW)', 'F_Xiph_WeightNormalized_100_Swing(%N/BW)', 'F_Umbi_WeightNormalized_100_Stride(%N/BW)', 'F_Umbi_WeightNormalized_100_Contact(%N/BW)', 'F_Umbi_WeightNormalized_100_Swing(%N/BW)', 'F_UL_WeightNormalized_100_Stride(%N/BW)', 'F_UL_WeightNormalized_100_Contact(%N/BW)', 'F_UL_WeightNormalized_100_Swing(%N/BW)', 'F_LL_WeightNormalized_100_Stride(%N/BW)', 'F_LL_WeightNormalized_100_Contact(%N/BW)', 'F_LL_WeightNormalized_100_Swing(%N/BW)', 'F_FT_WeightNormalized_100_Stride(%N/BW)', 'F_FT_WeightNormalized_100_Contact(%N/BW)', 'F_FT_WeightNormalized_100_Swing(%N/BW)'};

Export_Excel_Curves100_Kinetic_Normalized_Weight = [Header_Curves100_Kinetic_Normalized_Weight; num2cell(Data_Curves100_Kinetic_Normalized_Weight)];
xlswrite (Output_File_Path_Full_Curves100_Kinetic_Normalized_Weight, Export_Excel_Curves100_Kinetic_Normalized_Weight); 


% NORMALIZED CURVES BY SEGMENT FRONTAL AREA   
Output_File_Path_Full_Curves100_Kinetic_Normalized_Area = [Export_File_Path '\' File_Name '_CURVES100_KINETIC_AREA_NORMALIZED' '.xls'];
      
Data_Curves100_Kinetic_Normalized_Area = [F_Xiph_AreaNormalized_100_Stride F_Xiph_AreaNormalized_100_Contact F_Xiph_AreaNormalized_100_Swing F_Umbi_AreaNormalized_100_Stride F_Umbi_AreaNormalized_100_Contact F_Umbi_AreaNormalized_100_Swing F_UL_AreaNormalized_100_Stride F_UL_AreaNormalized_100_Contact F_UL_AreaNormalized_100_Swing F_LL_AreaNormalized_100_Stride F_LL_AreaNormalized_100_Contact F_LL_AreaNormalized_100_Swing F_FT_AreaNormalized_100_Stride F_FT_AreaNormalized_100_Contact F_FT_AreaNormalized_100_Swing];
Header_Curves100_Kinetic_Normalized_Area = {'F_Xiph_AreaNormalized_100_Stride(N/m2)', 'F_Xiph_AreaNormalized_100_Contact(N/m2)', 'F_Xiph_AreaNormalized_100_Swing(N/m2)', 'F_Umbi_AreaNormalized_100_Stride(N/m2)', 'F_Umbi_AreaNormalized_100_Contact(N/m2)', 'F_Umbi_AreaNormalized_100_Swing(N/m2)', 'F_UL_AreaNormalized_100_Stride(N/m2)', 'F_UL_AreaNormalized_100_Contact(N/m2)', 'F_UL_AreaNormalized_100_Swing(N/m2)', 'F_LL_AreaNormalized_100_Stride(N/m2)', 'F_LL_AreaNormalized_100_Contact(N/m2)', 'F_LL_AreaNormalized_100_Swing(N/m2)', 'F_FT_AreaNormalized_100_Stride(N/m2)', 'F_FT_AreaNormalized_100_Contact(N/m2)', 'F_FT_AreaNormalized_100_Swing(N/m2)'};

Export_Excel_Curves100_Kinetic_Normalized_Area = [Header_Curves100_Kinetic_Normalized_Area; num2cell(Data_Curves100_Kinetic_Normalized_Area)];
xlswrite (Output_File_Path_Full_Curves100_Kinetic_Normalized_Area, Export_Excel_Curves100_Kinetic_Normalized_Area); 


%% SECTION 16
% GRAPHICS AND GIFS %


    % SEGMENTS ANGULAR GRAPHICS
    
figure('Name','Segments Angles')
title('Angles');
hold on
plot(Foot_Angle(:,1),'r','LineWidth',2.5)
hold on
plot(Foot_Angle_100_Stride,'r:', 'LineWidth', 3.5)
hold on
plot (Shank_Angle(:,1),'b', 'LineWidth',2.5)
hold on
plot(Shank_Angle_100_Stride,'b:', 'LineWidth',3.5)
hold on
plot (Thigh_Angle(:,1), 'k','LineWidth',2.5)
hold on
plot (Thigh_Angle_100_Stride(:,1), 'k:','LineWidth',3.5)
hold on
plot (Trunk_Angle(:,1), 'g','LineWidth',2.5)
hold on
plot (Trunk_Angle_100_Stride(:,1), 'g:','LineWidth',3.5)
%plot([TO TO],[0 100],'-m', 'LineWidth',2.5)
hold on
%plot([0 120],[90 90],'-k', 'LineWidth',2.5)
legend ('Foot', 'Foot 100', 'Shank', 'Shank 100', 'Thigh', 'Thigh 100' , 'Trunk', 'Trunk 100');


figure('Name','Foot Segment Angles 100%')
title('Angles');
hold on
plot(Foot_Angle_100_Stride,'r:', 'LineWidth', 3.5)
hold on
plot(Foot_Angle_100_Contact,'b:', 'LineWidth', 3.5)
hold on
plot(Foot_Angle_100_Swing,'k:', 'LineWidth', 3.5)
legend ('Foot Stride 100', 'Foot Contact 100', 'Foot Swing 100');


figure('Name','Shank Segment Angles 100%')
title('Angles');
hold on
plot(Shank_Angle_100_Stride,'r:', 'LineWidth', 3.5)
hold on
plot(Shank_Angle_100_Contact,'b:', 'LineWidth', 3.5)
hold on
plot(Shank_Angle_100_Swing,'k:', 'LineWidth', 3.5)
legend ('Shank Stride 100', 'Shank Contact 100', 'Shank Swing 100');


figure('Name','Thigh Segment Angles 100%')
title('Angles');
hold on
plot(Thigh_Angle_100_Stride,'r:', 'LineWidth', 3.5)
hold on
plot(Thigh_Angle_100_Contact,'b:', 'LineWidth', 3.5)
hold on
plot(Thigh_Angle_100_Swing,'k:', 'LineWidth', 3.5)
legend ('Thigh Stride 100', 'Thigh Contact 100', 'Thigh Swing 100');


figure('Name','Shank-Foot Segment Angles 100%')
title('Angles')
xlabel('Shank Angular Position (º)');
ylabel ('Foot Angular Position (º)');
hold on
plot(Shank_Angle_100_Stride, Foot_Angle_100_Stride,'b:', 'LineWidth', 3.5)
hold on
plot(Shank_Angle_100_Stride(1), Foot_Angle_100_Stride(1),'ro', 'LineWidth', 3.5)
legend ('Shank-Foot Angles','TD')


figure('Name','Thigh-Shank Segment Angles 100%')
title('Angles')
xlabel('Thigh Angular Position (º)');
ylabel ('Shank Angular Position (º)');
hold on
plot(Thigh_Angle_100_Stride, Shank_Angle_100_Stride,'b:', 'LineWidth', 3.5)
hold on
plot(Thigh_Angle_100_Stride(1), Shank_Angle_100_Stride(1),'ro', 'LineWidth', 3.5)
legend ('Thigh-Shank Angles','TD')


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
plot(Shank_Angle_100_Stride_Normalized(1:TO_in_Percentage), Shank_Angular_Speed_100_Stride_Normalized(1:TO_in_Percentage), 'g','LineWidth',2.5)
hold on
plot(Shank_Angle_100_Stride_Normalized(TO_in_Percentage:end), Shank_Angular_Speed_100_Stride_Normalized(TO_in_Percentage:end), 'm','LineWidth',2.5)
hold on
plot(Shank_Angle_100_Stride_Normalized(1),Shank_Angular_Speed_100_Stride_Normalized(1), 'k*', 'LineWidth', 3.0)
hold on
plot(Thigh_Angle_100_Stride_Normalized(1:TO_in_Percentage), Thigh_Angular_Speed_100_Stride_Normalized(1:TO_in_Percentage), 'b','LineWidth',2.5)
hold on
plot(Thigh_Angle_100_Stride_Normalized(TO_in_Percentage:end), Thigh_Angular_Speed_100_Stride_Normalized(TO_in_Percentage:end), 'r','LineWidth',2.5)
hold on
plot(Thigh_Angle_100_Stride_Normalized(1),Thigh_Angular_Speed_100_Stride_Normalized(1), 'k*', 'LineWidth', 3.0)
legend ('Shank Contact','Shank Swing', 'TD', 'Thigh Contact', 'Thigh Swing', 'TD')


figure('Name','Phase Angle and Relative Phase Foot/Shank')
title('Phase Angle');
xlabel('Frame (n)');
ylabel ('Phase Angle (º)');
hold on
plot(Foot_Phase_Angle(:,1), 'r','LineWidth',2.5)
hold on
plot(Shank_Phase_Angle(:,1), 'g','LineWidth',2.5)
hold on
plot (Relative_Phase_Foot_Shank(:,1), 'b', 'LineWidth', 2.5)
legend('Foot','Shank', 'Relative Phase Foot/Shank')


figure('Name','Phase Angle and Relative Phase Shank/Thigh')
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



figure('Name','Vector Coding Shank Foot')
title('Vector Coding Shank Foot STRIDE');
xlabel('Frame (n)');
ylabel ('Coupling Angle (º)');
plot(Coupling_Angle_Shank_Foot_Stride, 'g','LineWidth',2.5)
legend('Coupling Anlge Shank Foot')


figure('Name','Vector Coding Thigh Shank')
title('Vector Coding Thigh Shank STRIDE');
xlabel('Frame (n)');
ylabel ('Coupling Angle (º)');
plot(Coupling_Angle_Thigh_Shank_Stride, 'g','LineWidth',2.5)
legend('Coupling Anlge Thigh Shank')


        % JOINTS ANGULAR GRAPHICS
           
figure('Name','Ankle Joint Angles 100%')
title('Angles');
hold on
plot(Ankle_Joint_Angle_100_Stride,'r:', 'LineWidth', 3.5)
hold on
plot(Ankle_Joint_Angle_100_Contact,'b:', 'LineWidth', 3.5)
hold on
plot(Ankle_Joint_Angle_100_Swing,'k:', 'LineWidth', 3.5)
legend ('Ankle Joint Stride 100', 'Ankle Joint Contact 100', 'Ankle Joint Swing 100');


figure('Name','Knee Joint Angles 100%')
title('Angles');
hold on
plot(Knee_Joint_Angle_100_Stride,'r:', 'LineWidth', 3.5)
hold on
plot(Knee_Joint_Angle_100_Contact,'b:', 'LineWidth', 3.5)
hold on
plot(Knee_Joint_Angle_100_Swing,'k:', 'LineWidth', 3.5)
legend ('Knee Joint Stride 100', 'Knee Joint Contact 100', 'Knee Joint Swing 100');


figure('Name','Hip Joint Angles 100%')
title('Angles');
hold on
plot(Hip_Joint_Angle_100_Stride,'r:', 'LineWidth', 3.5)
hold on
plot(Hip_Joint_Angle_100_Contact,'b:', 'LineWidth', 3.5)
hold on
plot(Hip_Joint_Angle_100_Swing,'k:', 'LineWidth', 3.5)
legend ('Hip Joint Stride 100', 'Hip Joint Contact 100', 'Hip Joint Swing 100');


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



        % DRAG FORCE GRAPHICS

figure ('Name', 'Drag Force Trunk')
title ('Drag Force Trunk');
xlabel('Frames (n)');
ylabel ('Drag Force (N)');
hold on
plot(F_Xiph,'r','LineWidth',1.5);
hold on
plot(F_Xiph_100_Stride,'r:','LineWidth',2.0);
hold on
plot(F_Xiph_100_Contact,'b','LineWidth',2.0);
hold on;
plot(F_Xiph_100_Swing,'k','LineWidth',2.0);
legend ('F_Xiph', 'F_Xiph_100_Stride','F_Xiph_100_Contact','F_Xiph_100_Swing');

figure ('Name', 'Drag Force Upper Leg and Torque')
title ('Drag Force Upper Leg');
xlabel('Frames (n)');
ylabel ('Drag Force (N)');
hold on
plot(F_UL,'r','LineWidth',1.5);
hold on
plot(F_UL_100_Stride,'r:','LineWidth',2.0);
hold on
plot(F_UL_100_Contact,'b','LineWidth',2.0);
hold on;
plot(F_UL_100_Swing,'k','LineWidth',2.0);
hold on
plot(T_UL, 'go', 'LineWidth',3.5)
legend ('F_UL', 'F_UL_100_Stride','F_UL_100_Contact','F_UL_100_Swing','UL Torque');


figure ('Name', 'Drag Force Lower Leg and Torque')
title ('Drag Force Lower Leg');
xlabel('Frames (n)');
ylabel ('Drag Force (N)');
hold on
plot(F_LL,'r','LineWidth',1.5);
hold on
plot(F_LL_100_Stride,'r:','LineWidth',2.0);
hold on
plot(F_LL_100_Contact,'b','LineWidth',2.0);
hold on;
plot(F_LL_100_Swing,'k','LineWidth',2.0);
hold on
plot(T_LL, 'go', 'LineWidth',3.5)
legend ('F_LL', 'F_LL_100_Stride','F_LL_100_Contact','F_LL_100_Swing', 'LL Torque');


figure ('Name', 'Drag Force Foot and Torque')
title ('Drag Force Foot');
xlabel('Frames (n)');
ylabel ('Drag Force (N)');
hold on
plot(F_FT,'r','LineWidth',1.5);
hold on;
plot(F_FT_100_Stride,'r:','LineWidth',2.0);
hold on;
plot(F_FT_100_Contact,'b','LineWidth',2.0);
hold on;
plot(F_FT_100_Swing,'k','LineWidth',2.0);
hold on
plot(T_FT, 'go', 'LineWidth',3.5)
legend ('F_FT', 'F_FT_100_Stride','F_FT_100_Contact','F_FT_100_Swing', 'FT Torque');


    % GIF WITH GAIT VIDEO CREATION

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