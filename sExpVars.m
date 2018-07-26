%  'RestrictProcessing' Restrict stimulus processing to a specific subarea
%     of the screen. If your visual stimulus only covers a subarea of the
%     display screen you can restrict PTB's output processing to that
%% experimental variables
%%%%%%%%%% Exp. vars
trl.TrainingTrials=5;
trl.ntrl = 100;
trl.num = 0;

%% rand number stream init
if verLessThan('matlab', '8.0.0')
    rand('seed', sum(100*clock)); 
else    
    RandStream('mt19937ar', 'Seed', sum(100*clock));        
end

%% Init graphical variables (what defines stimulus)
[gfx.scMid(1) gfx.scMid(2)] = RectCenter(winRect);
gfx.ifi=Screen('GetFlipInterval', win); % in seconds

gfx.screen = 'Dinamondtron NF, Gatweay Vx920 17in';

gfx.monitor_freq    = 85; % Hz
gfx.ms_per_frame	= 1000/gfx.monitor_freq; 
gfx.h_monitor_size  = 36; % cm
gfx.v_monitor_size  = 27.5; % cm

gfx.h_pixel         = winRect(3); 
gfx.v_pixel         = winRect(4); 

%gfx.scMid(1) = gfx.h_pixel*.5;
%gfx.scMid(2) = gfx.v_pixel*.5;

gfx.sub_distance    = 57; % distance of subject from screen cm
gfx.h_deg2pix       = 1/((180/pi) * atan((gfx.h_monitor_size/gfx.h_pixel) / gfx.sub_distance));
gfx.v_deg2pix       = 1/((180/pi) * atan((gfx.v_monitor_size/gfx.v_pixel) / gfx.sub_distance));
    

gfx.cohere = 1;

% Colormap
gfx.black = BlackIndex(win); 
gfx.white = WhiteIndex(win);
gfx.gray =  round((gfx.black + gfx.white) / 2);      
gfx.inc = gfx.white - gfx.gray;

util.direxp = cd;

% Dots parameters
gfx.dir  = 0;

gfx.h_pix2deg = 1/gfx.h_deg2pix;
gfx.v_pix2deg = 1/gfx.v_deg2pix;

gfx.fixsize = 12;                                                          % Fixation dot size

gfx.red=[255 0 0];
gfx.blue=[0 0 255];
gfx.green=[0 255 0];

gfx.height = 4;                                                            
gfx.width = 30;                                                            
gfx.height_pix = gfx.height./gfx.v_pix2deg;
gfx.width_pix = gfx.width./gfx.h_pix2deg;

% dial display
gfx.Let{4} = 'qwertyuiop';                                             % QWERTY 
gfx.Let{3} = 'asdfghjkl';                            
gfx.Let{2} = 'zxcvbnm<';                              
gfx.Let{1} = '___';                                                    % spacebar ... 

% Magnification factor
gfx.M = 1.0;                                                        

ly=1;                                                               
gfx.bE(ly) = fix(20*gfx.h_deg2pix);                                    % Lower eccentricity limit of bounding box
gfx.bD(ly) = fix(2.2*gfx.h_deg2pix);                                   % Box radial length
gfx.ArcDeg = 5*(pi/180);                                               % Arc width (deg) for each box
gfx.bDa(ly) = gfx.ArcDeg;

nletters=length(gfx.Let{1});                                           % angles depending on the ly
gfx.maxAngle(1) = nletters.*gfx.bDa(1);
gfx.b(1).C = linspace(0,gfx.maxAngle(1), nletters)-gfx.maxAngle(1)*.5; % angular centre of the boxes

% for next lys, E calculated from D
for ly = 2:4
    nletters=length(gfx.Let{ly});
    gfx.bE(ly) = (gfx.bD(ly-1)+gfx.bE(ly-1)).*gfx.M; 
    gfx.bD(ly) = gfx.bE(ly)-gfx.bE(ly-1);   
    gfx.bDa(ly) = gfx.ArcDeg.*(gfx.bD(ly)/gfx.bD(ly-1)); % arc size
    gfx.maxAngle(ly) = nletters.*gfx.bDa(ly);
    
    gfx.b(ly).C = linspace(0,gfx.maxAngle(ly), nletters)-gfx.maxAngle(ly)*.5;        % angular centre of the boxes
end
ly = 5;
gfx.bE(ly) = (gfx.bD(ly-1)+gfx.bE(ly-1)).*gfx.M; 
gfx.b(ly).C = gfx.b(ly-1).C; 

% CENTRE OF THE KEYBOARD CIRCLE
gfx.kbCentre = [gfx.scMid(1) gfx.scMid(2)+750];
gfx.wrdCentre = [gfx.scMid(1) gfx.scMid(2)+250]; %gfx.Box(ly).E(1)];

%% log file name
util.dir_name = sprintf('%s/data/e%ds%d/',util.direxp, trl.exp, trl.sub);
if ~exist(util.dir_name,'dir')
    mkdir(util.dir_name);
end
util.log_name = sprintf('%se%ds%d.log',util.dir_name, trl.exp, trl.sub);
