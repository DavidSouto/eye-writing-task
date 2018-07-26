    %% initialize header
    Disp.name = strcat('File: ', Man.currfile, ', Sub: ', int2str(DATA(Trl.trl,Col.sub)), ', Err: ', int2str(error(Trl.num)));

    %% vector with time, scale for axis
    Disp.TIME   = Ctrl.sdur: Ctrl.sdur : length(X);
    TEMP        = Disp.TIME';
    DIR = [1 -1];

    %% setup plot
    figure(Man.plot);
    cla reset;

%    subplot(311);
    cla reset;
    title(['Velocity:' Disp.name]);

    crit_interval = [stimOn 700];
 
    dirid = DATA(Trl.num,Col.sideidx);   
    hold on; plot([TEMP(stimOn) TEMP(stimOff)], Trl.SIDE(dirid+1).*[Ctrl.amplitude Ctrl.amplitude],'b-');
    hold on; plot([TEMP(stimOn) TEMP(stimOff)], [0 0],'--k');
 
    display_int = stimOn:stimOff;
    hold on; plot(TEMP(display_int),Y(display_int),'b'); % not dDY
    hold on; plot(TEMP(display_int),X(display_int),'r'); % not dDY
    
    % plot saccade on- and offsets
    if ~isempty(SACCADE)
        hold on; plot(Ctrl.sdur * SACCADE(:,1), X(SACCADE(:,1)), '+r');
        hold on; plot(Ctrl.sdur * SACCADE(:,2), X(SACCADE(:,2)), '+g');
        
        hold on; plot(Ctrl.sdur * SACCADE(:,1), Y(SACCADE(:,1)), '+r');
        hold on; plot(Ctrl.sdur * SACCADE(:,2), Y(SACCADE(:,2)), '+g');
    end
    
    text(TEMP(stimOff)-500,Man.PosRange(1)-1,sprintf('Instruction: %s',Trl.ANT{DATA(Trl.num,Col.antidx)+1}));
    
    % plot bars and trial info
    sPlotEvents(Ctrl, EVE, Man.off, Man.PosRange);
    Trl.xstr    = stimOn;
    Trl.ystr    = 0;
    sPlotPar(Trl, Col, Err, error(Trl.num));

    ylabel('Position'); xlabel('Time (ms)');
    axis([stimOn stimOff -Man.PosRange Man.PosRange]);
    %set(gca, 'XTick', 0:200:3000, 'XTickLabel', 0:200:3000-stimOn); 
    