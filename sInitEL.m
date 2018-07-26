%% Init eyelink variables
el = EyelinkInitDefaults(win); 
el.backgroundcolour = 0;
el.foregroundcolour = 128;   
el.LOSTDATAEVENT = hex2dec('3F');
el.msgfont = 'Arial';
el.msgfontsize = 40;
el.msgfontcolour = [0 0 0];
el.backgroundcolour = 128;
el.foregroundcolour = 255;   
el.calibrationtargetsize = 1;
el.calibrationtargetwidth = .1;
el.calibrationtargetcolour = 0;

EyelinkUpdateDefaults(el);

if ~EyelinkInit(0, 1)
     fprintf('Eyelink Init aborted.\n');
     cleanup;  % cleanup function
     return;
end

%util.el_name   = ['s',int2str(trl.sub),'t',int2str(1),'.edf'];  
util.el_name   = ['s',int2str(trl.sub),'t',int2str(trl.blk),'.edf'];                

elopen = Eyelink('Openfile', util.el_name);
if elopen
    error('edf file failed to open');
    cleanup;
end                    

Eyelink('command', 'active_eye = RIGHT');

%% DEFINE CALIBRATION
Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, gfx.scMid(1)*2-1, gfx.scMid(2)*2-1);
Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, gfx.scMid(1)*2-1, gfx.scMid(2)*2-1);
    
trl.full_calib = 1;

if trl.full_calib==0
        
    nsamples = 9;

    x = gfx.scMid(1);
    y = gfx.scMid(2);
    grid_width = 8;
    grid_height = 5;
    
    x_off = round(grid_width*gfx.h_deg2pix);
    y_off = round(grid_height*gfx.v_deg2pix);
    calib = sprintf('%d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d',...
        x,y,...
        x,y-y_off,...
        x,y+y_off,...
        x-x_off,y,...
        x+x_off,y,...
        x-x_off,y-y_off,...
        x+x_off,y-y_off,...
        x-x_off,y+y_off,...
        x+x_off,y+y_off);
       
elseif trl.full_calib == 1
    
    nsamples = 13;

    % check pedro's custom calibration
    
    % [width, height]=Screen('WindowSize', screenNumber);
    %Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, width-1, height-1);
    %Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, width-1, height-1);
    x = gfx.scMid(1);
    y = gfx.scMid(2);
    grid_hwidth = 4; %half width
    grid_hheight = 2.5; %half height
    
    x_off = round(grid_hwidth*gfx.h_deg2pix);
    y_off = round(grid_hheight*gfx.v_deg2pix);
    
    % 13-point calibration, meaning 1 in the centre and a 4 by 3 grid         
    %[X,Y]=meshgrid(linspace(x-x_off,x+x_off,4),linspace(y-y_off,y+y_off,3));
    X(1,1:4) = linspace(x-x_off,x+x_off,4);  % 4, 5, 4
    X(2,1:5) = linspace(x-x_off,x+x_off,5);
    X(3,1:4) = linspace(x-x_off,x+x_off,4);

    Y(1,:) = linspace(y-y_off,y+y_off,3);  % 4, 5, 4

    calib = sprintf('%d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d, %d %d, %d %d, %d %d, %d %d, %d',...
        X(1,1),Y(1),... % First row
        X(1,2),Y(1),... %
        X(1,3),Y(1),... %
        X(1,4),Y(1),... %
        X(2,1),Y(2),... % Second row
        X(2,2),Y(2),... %
        X(2,3),Y(2),... %
        X(2,4),Y(2),... %
        X(2,5),Y(2),... %
        X(2,1),Y(3),... % Third row
        X(2,2),Y(3),... %
        X(2,3),Y(3),... %
        X(2,4),Y(3));%,... %
    
end

Eyelink('command',['calibration_samples = ',int2str(nsamples)]);
Eyelink('command',['calibration_type = HV' int2str(nsamples)]);
Eyelink('command','generate_default_targets = NO');

Eyelink('command',sprintf('calibration_targets = %s',calib));
Eyelink('command',sprintf('validation_targets = %s',calib));
Eyelink('command','validation_sequence = 0,1,2,3,4,5,6,7,8,9,10,11,12,13');
Eyelink('command','calibration_sequence = 0,1,2,3,4,5,6,7,8,9,10,11,12,13');
Eyelink('command','button_function 1 ''accept_target_fixation''');    

% set parser: psychophysics recommended settings (high sensitivity)
Eyelink('command', 'recording_parse_type = GAZE');
Eyelink('command', 'saccade_velocity_threshold = 22');
Eyelink('command', 'saccade_acceleration_threshold = 5000');
Eyelink('command', 'saccade_motion_threshold = 0.0');
Eyelink('command', 'saccade_pursuit_fixup = 60');
Eyelink('command', 'fixation_update_interval = 0');
Eyelink('command', 'set_cal_sounds','off');
Eyelink('command', 'set_dcorr_sounds','off');

% use gamepad to accept fixation or do drift-correction
Eyelink('command', 'button_function 5 "accept_target_fixation"');

[result, reply]=Eyelink('ReadFromTracker','enable_automatic_calibration');
    
if reply % reply = 1
        fprintf('Automatic sequencing ON');
else
        fprintf('Automatic sequencing OFF');
end


% Can't use this with ELII

Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,FIXATION,SACCADE,AREA,GAZERES,HREF,PUPIL,STATUS,INPUT,HMARKER')


% set EDF file contents
% Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON'); 
% Eyelink('command', 'file_sample_data  = LEFT, RIGHT, GAZE, AREA, SACCADE, BLINK, MESSAGE'); 
Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON'); 
Eyelink('command', 'file_sample_data  = LEFT, RIGHT, GAZE, AREA, SACCADE, BLINK, MESSAGE'); 

%if Eyelink('command','inputword_is_window = ON')
%    error('inputword_is_window error')
%end

% FIX SAMPLING RATE!
Eyelink('command', 'sample_rate = %d',500);
