if trl.Eyelink,
    Screen('Flip',win);
    Eyelink('StopRecording');                
    % record data to txt file
    status = EL2_playback(util.el_datdir);
    if ~status;
        error('playback failed');
    end
end