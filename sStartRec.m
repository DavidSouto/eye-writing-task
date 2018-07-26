% Open eyelink data file and check recording
if trl.Eyelink
    % [util] = sStartRec;
%    util.el_datname   = ['e', int2str(trl.exp),'s',int2str(trl.sub),'b',int2str(trl.blk),'t',int2str(trl.num),'.dat'];
%    util.el_datdir = [cd '\data\e', int2str(trl.exp),'s',int2str(trl.sub),'\' util.el_datname]; 

%    Eyelink('StartRecording');

    % [util] = sStartRec;
    util.el_datname   = ['b',int2str(trl.blk),'t',int2str(trl.num),'.edf'];

    if exist([util.dir_name util.el_datname],'file'),

        sca;
        error('THIS EDF FILE BLOCK NUMBER ALREADY EXISTS');
        
    end
    
    elopen = Eyelink('Openfile', [util.dir_name util.el_datname]);
    Eyelink('command', 'active_eye = RIGHT');

    if elopen,
        error('edf file failed to open');
        cleanup;
    end                    

    Eyelink('StartRecording');    
    
    el.eye_used = Eyelink('EyeAvailable');
    switch el.eye_used
        case el.BINOCULAR
            error('tracker indicates binocular');
        case el.LEFT_EYE
            error('tracker indicates left eye');
        case el.RIGHT_EYE
            disp('tracker indicates right eye');
        case -1
            error('eye available returned -1');
        otherwise
            error('unexpected result from eye available');
    end
end          