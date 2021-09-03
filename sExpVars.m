%  'RestrictProcessing' Restrict stimulus processing to a specific subarea
%     of the screen. If your visual stimulus only covers a subarea of the
%     display screen you can restrict PTB's output processing to that
%% experimental variables
%%%%%%%%%% Exp. vars
trl.TrainingTrials=5;
trl.ntrl = 100;
trl.num = 1;

%% rand number stream init
if verLessThan('matlab', '8.0.0')
    rand('seed', sum(100*clock)); 
else    
    RandStream('mt19937ar', 'Seed', sum(100*clock));        
end

%% Init graphical variables (what defines stimulus)
[gfx.scMid(1), gfx.scMid(2)] = RectCenter(winRect);
gfx.ifi=Screen('GetFlipInterval', win); % in seconds

gfx.screen = 'GDC 5.05';

gfx.monitor_freq    = 75; % Hz
gfx.ms_per_frame	= 1000/gfx.monitor_freq; 

gfx.h_monitor_size  = 40; % cm
gfx.v_monitor_size  = 30.0; % cm

gfx.h_pixel         = 1280; % 1028
gfx.v_pixel         = 1024; % 768
gfx.sub_distance    = 61;   % distance of subject from screen cm
gfx.h_deg2pix       = 1/((180/pi) * atan((gfx.h_monitor_size/gfx.h_pixel) / gfx.sub_distance));
gfx.v_deg2pix       = 1/((180/pi) * atan((gfx.v_monitor_size/gfx.v_pixel) / gfx.sub_distance));

gfx.cohere = 1;

% Colormap
gfx.black = BlackIndex(win); 
gfx.white = WhiteIndex(win);
gfx.gray =  round((gfx.black + gfx.white) / 2);      
gfx.inc = gfx.white - gfx.gray;

util.direxp = '..';

gfx.cohere = 1; %.5;

% Dots parameters
gfx.dir  = 0;

gfx.h_pix2deg = 1/gfx.h_deg2pix;
gfx.v_pix2deg = 1/gfx.v_deg2pix;

gfx.fixsize = 12; % Fixation dot size

gfx.red=[255 0 0];
gfx.blue=[0 0 255];
gfx.green=[0 255 0];

gfx.height = 4;                                                            
gfx.width = 30;                                                            
gfx.height_pix = gfx.height./gfx.v_pix2deg;
gfx.width_pix = gfx.width./gfx.h_pix2deg;

gfx.WordDur = .5; % display the target word for 

trl.NWORDS = 100;

% QWERTY
gfx.kLetter{3} = 'qwertyuiop';                                        
gfx.kLetter{2} = 'asdfghjkl';                            
gfx.kLetter{1} = 'zxcvbnm';    

% square size
gfx.vkeySize = 95;

gfx.ncols = 10; % columns in the vkeyboard
gfx.nrows = 3; % rows in the vkeyboard

% Position of the vkeyboard lines, centered on the screen centre
hhalf = .5*gfx.ncols*gfx.vkeySize;
gfx.x_vkeyLine = linspace(-hhalf, hhalf, gfx.ncols+1);

vhalf = .5*gfx.nrows*gfx.vkeySize;
gfx.y_vkeyLine = linspace(-vhalf, vhalf, gfx.nrows+1);

% Get the center of the squares, where letters are positioned
gfx.x_vkeyCentre = gfx.x_vkeyLine+gfx.vkeySize*.5;
gfx.y_vkeyCentre = gfx.y_vkeyLine+gfx.vkeySize*.5;

% QWERTY vkeyboard
gfx.Letter{1} = 'qwertyuiop';                                               
gfx.Letter{2} = 'asdfghjkl ';                            
gfx.Letter{3} = ' zxcvbnm  ';

% centering of words to be written
gfx.wrdCentre = [gfx.scMid(1) gfx.scMid(2)+200]; 
        
%% log file name
util.dir_name = sprintf('%s/data/e%ds%d',util.direxp, trl.exp, trl.sub);
if ~exist(util.dir_name,'dir')
    mkdir(util.dir_name);
end
util.log_name = sprintf('%s/e-e%ds%d.log',util.dir_name, trl.exp, trl.sub);
