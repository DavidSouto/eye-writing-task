function fCloseEL(util)
    Eyelink('command', 'generate_default_targets = YES');
    Eyelink('CloseFile');
    Eyelink('ReceiveFile',[],util.dir_name,1);  
    Eyelink('Shutdown');
    sca;
    commandwindow; 
end
