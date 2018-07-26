if trl.Eyelink,
    %Screen('Flip',win);
    Eyelink('StopRecording');                

    % record data to txt file (would need a .mexa64 version)
    %status = EL2_playback(util.el_datdir); % did not compile under linux 64bit   
    %if ~status;
    %    error('playback failed');
    %end
    
    Eyelink('CloseFile');
    Eyelink('ReceiveFile',[],util.dir_name,1);  
    
end