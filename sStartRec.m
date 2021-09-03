% Open eyelink data file and check recording
if trl.Eyelink,
    % [util] = sStartRec;
    util.el_datname   = ['e-b',int2str(trl.blk),'t',int2str(trl.num),'.dat'];
    util.el_datdir = [pwd '\..\data\e', int2str(trl.exp),'s',int2str(trl.sub),'\e-b',int2str(trl.blk),'t',int2str(trl.num),'.dat']; 
    Eyelink('StartRecording');

    el.eye_used = Eyelink('EyeAvailable');
    switch el.eye_used,
        case el.BINOCULAR
            error('tracker indicates binocular');
        case el.LEFT_EYE
            disp('tracker indicates left eye');
        case el.RIGHT_EYE
            error('tracker indicates right eye');
        case -1
            error('eye available returned -1');
        otherwise
            error('unexpected result from eye available');
    end
end          