% need to start by taking an array of x,y raw eye movements/mouse movements
% and test an algorithm to classify fixations on-line that works as well for both

% Try algorithm downloaded from
% http://www.humlab.lu.se/en/person/MarcusNystrom/ Nystrom and Holmqvist
% adaptive algorithm
addpath(genpath('./FixationDetector'))

% Could have a baseline period, during which fixation detection thresholds
% are set ...

x = rand(100,1);
y = rand(100,1);

% Need to give a try at this:
% I-VS with Kalman filter wtih realistic eyeball characteristics; shown to work with EL1000: http://cs.txstate.edu/~ok11/igaze_emd_online.html
% iGaze.rar

% Could use own: Vel + Cluster + Dur + vel? Points close together ... 
