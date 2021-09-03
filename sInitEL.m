%% Init eyelink variables
el = EyelinkInitDefaults(win);   
el.LOSTDATAEVENT = hex2dec('3F');
el.msgfont = 'Arial';
el.msgfontsize = 40;
el.msgfontcolour = [0 0 0];
el.backgroundcolour = 128;
el.foregroundcolour = 0;   
el.calibrationtargetsize = 1.5;
el.calibrationtargetwidth = .4;
el.calibrationtargetcolour = 0;

EyelinkUpdateDefaults(el);

if ~EyelinkInit(0, 1)
     fprintf('Eyelink Init aborted.\n');
     cleanup;  % cleanup function
     return;
end

%util.el_name   = ['s',int2str(trl.sub),'t',int2str(1),'.edf'];  
util.el_name   = ['s',int2str(trl.sub),'b',int2str(trl.blk),'.edf'];                

elopen = Eyelink('Openfile', util.el_name);
if elopen
    error('edf file failed to open');
    cleanup;
end                    

Eyelink('command', 'active_eye = LEFT');

%% DEFINE CALIBRATION
x = gfx.scMid(1);
y = gfx.scMid(2);
x_off = round(gfx.CALIB_X*gfx.h_deg2pix);
y_off = round(gfx.CALIB_Y*gfx.v_deg2pix);
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

Eyelink('command','calibration_type = HV9');
Eyelink('command','generate_default_targets = NO');
Eyelink('command',sprintf('calibration_targets = %s',calib));
Eyelink('command',sprintf('validation_targets = %s',calib));
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

Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,FIXATION,SACCADE,AREA,GAZERES,HREF,PUPIL,STATUS,INPUT,HMARKER')

Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON'); 
Eyelink('command', 'file_sample_data  = LEFT, RIGHT, GAZE, AREA, SACCADE, BLINK, MESSAGE'); 

if Eyelink('command','inputword_is_window = ON')
    error('inputword_is_window error')
end

% FIX SAMPLING RATE!
Eyelink('command', 'sample_rate = %d',1000);
