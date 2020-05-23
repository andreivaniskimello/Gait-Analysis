# Gait-Analysis
This is MatLab Routine to Gait Analysis.
It have been developted to Shallow Water Walking.

There are 04 routines. Each one is for one specific water depth: Knee, Hip, Umbilical, Xiphoid.
This distinction is made to account the different immersed body part in each depth, in order to adjust the proper Drag Force Model.

The routines are made to run with kinematic data from 01 immersed camera in the sagital plane. 

INPUTS are:   
    % POSITION x TIME Matrix of LOWER LIMBR AND TRUNK POINTS
    % SUBJECT BODY MASS
    % IMMERSION DEPTH
    % Touch-Down and Take-Off Frames of ONE STRIDE
    % LOW-PASS BUTTERWORTH FILTER PARAMETERS
    % SAMPLE FREQUENCY
    % NUMBER OF INTERPOLATION POINTS TO NORMALIZE CURVES TO 0 - 100%
    % STRIDE, CONTACT, SWING PHASES DURATIONS
    
    OUTPUTS are:
    % SPATIOTEMPORAL DATA
    % ANGULAR (SEGMENTS AND JOINTS) DATA
    % INTRALIMB COORDINATION PARAMETERS
    % DRAG FORCE
    % GROUND REACTION FORCES
