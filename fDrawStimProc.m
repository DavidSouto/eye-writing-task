function [trl, POS]= fDrawStimProc(win,gfx,trl,util,el)
   
    if trl.Eyelink
        Eyelink('Message', 'TRIALID %d', trl.trl); % EVENT 1

        % This supplies the title at the bottom of the eyetracker display
        Eyelink('command', 'record_status_message "TRIAL %d/%d"', trl.trl, 100);%, char(gfx.imLUT.imname(IMTYPE+1))); 
    end
    WaitSecs(0.1);
    
    global vbl;
    vbl = Screen('Flip', win);

    global POS;                                                            % Global record of cursor position
    POS.x = NaN*ones(5000,1);                                              % max ~58 second buffer
    POS.y = NaN*ones(5000,1);
    POS.t = NaN*ones(5000,1);
    POS.EVT = NaN*ones(5000,1); % get character with char()
    
    global dwelling;
    global layerid;
    global letterid;
    global oletter;
    global ox; % old x and y
    global oy;
    ox = NaN;
    oy = NaN;
    dwelling = 0;                                                          % Color of cursor ...    word
    letterid = [];
    oletter = [];                                                          % old and new letters
    layerid = [];
    trl.ccolor = [0 0 0];                                                  % Number of samples dwelling on the 
    trl.alpha = 1; % proportion updated
    if trl.Eyelink                                                         % real-time mode
        Eyelink('command','begin_realtime_mode(100)');              
        Eyelink('Message','SYNCTIME'); % EVENT 2 / not always event 2 ...
    end
    
    priorityLevel=MaxPriority(win);
    Priority(priorityLevel);

	trl.complete = 0;
    trl.repeat = 0;
    trl.word = '';
    SetMouse(gfx.wrdCentre(1),gfx.wrdCentre(2));                           % Init mouse location at center
    
    try
        nDisplayKeyboard(win,gfx);                                         % dial virtual keyboard
        %Screen('DrawingFinished');
        %img = Screen('GetImage',win);
        %imwrite(img,'img.jpg')
    catch me
        sca
        keyboard
    end
    trl = nDisplayWord(win,gfx,trl);     % display word to write
    
    % starting to record the cursor from there
    if trl.joystickControl
        fWaitForJoyStick()
    else
        KbWait([], 2);
    end
    
    EVT = '!'; % trial start event
    t0 = GetSecs();
    trl=nGetCursorCoords(trl);                               % cursor x,y coordinates, from eye or mouse    
    trl=nTestSelection(gfx,trl,EVT);    
    while (GetSecs()-t0)<0.5 % this is to avoid the eyetracker gets stuck and we can't move on
        EVT = '0';
        trl=nTestSelection(gfx,trl,EVT);
        trl=nGetCursorCoords(trl);                           % cursor x,y coordinates, from eye or mouse
        nDrawFeedback(win,gfx,trl);                          % show cursor position
        if ~trl.complete
            t0 = GetSecs();
            trl= nGetResp(trl,el,util);                                             % opportunity to skip if incomplete
        end  
    end
    
    %criterion_filename = [int2str(trl.sub) '.savecrit.mat'];% Dwelling criterion 
    VAL = trl.DwellingCriterion;
    save(trl.criterion_filename,'VAL');
    
    fprintf('trl.crit %d\n', trl.DwellingCriterion);

    fprintf('trial %d repeat %d\n', trl.trl, trl.repeat);
    if trl.Eyelink
        Eyelink('Message','TRIAL_END');  % EVENT 3
        Eyelink('command','end_realtime_mode()'); 
    end
    Priority(0);
    
    % nested functions 
    function trl = nGetCursorCoords(trl)
        if trl.control          
            if trl.Eyelink
                % eye-control mode
                trl = nDetectFix(trl,el);  
            end
        else                                                               % mouse-control mode
            [x, y, ~] = GetMouse; 
            trl.x = x;            
            trl.y = y;
        end
                
        %% --- Exponential filtering ---
        a = 0.3;
        if isnan(ox)
            trl.fx = trl.x;
            trl.fy = trl.y;
            ox = trl.fx;
            oy = trl.fy;
        else
            trl.fx = a * trl.x + (1-a)*ox;
            trl.fy = a * trl.y + (1-a)*oy;    
            ox = trl.fx;
            oy = trl.fy;
        end

        %% --- Update position buffer with last sample ---
        POS.x(1:end-1) = POS.x(2:end);                                     
        POS.y(1:end-1) = POS.y(2:end);   
        POS.t(1:end-1) = POS.t(2:end);
        POS.x(end) = trl.x;
        POS.y(end) = trl.y;
        POS.t(end) = GetSecs();
    end
    function nDisplayKeyboard(win,gfx)
        
        %-------- Display the vkeyboard borders at screen centre (gfx.scMid)
        gfx.WordPixSize = 30;
        Screen('TextSize', win, gfx.WordPixSize);
        for row = 1:gfx.nrows+1 % horizontals
           Screen('DrawLines', win, [[gfx.x_vkeyLine(1);gfx.y_vkeyLine(row)], [gfx.x_vkeyLine(end);gfx.y_vkeyLine(row)]], 3, 0, gfx.scMid,2);        
        end
        for col = 1:gfx.ncols+1 % verticals

           Screen('DrawLines', win, [[gfx.x_vkeyLine(col);gfx.y_vkeyLine(1)], [gfx.x_vkeyLine(col);gfx.y_vkeyLine(end)]], 3, 0, gfx.scMid,2);        
        end

        %-------- Display the letters                                      
        for row = 1:gfx.nrows
            for col = 1:gfx.ncols

                 DrawFormattedText(win, gfx.Letter{row}(col), 'center', 'center', ...
                    [255 255 255],[],[],[],[],[],...
                    CenterRectOnPoint([0 0 gfx.WordPixSize gfx.WordPixSize], gfx.scMid(1)+gfx.x_vkeyCentre(col), gfx.scMid(2)+gfx.y_vkeyCentre(row)));   
            end
        end
        vbl = Screen('Flip', win, vbl + 0.5 * gfx.ifi, 1);                 
    end
    function trl = nTestSelection(gfx,trl,EVT)
        
       cols = [];
       rows = [];
       x= trl.fx - gfx.scMid(1);
       y= trl.fy - gfx.scMid(2);

       if(x < gfx.x_vkeyLine(end)) && (x > gfx.x_vkeyLine(1))
           cols = sum(x>gfx.x_vkeyLine);
       end
       if(y < gfx.y_vkeyLine(end)) && (y > gfx.y_vkeyLine(1))
           rows = sum(y > gfx.y_vkeyLine);
       end   
       if ~isempty(rows) && ~isempty(cols)                                 % Dwell time criterion
           trl.letter = gfx.Letter{rows}(cols);

           if  isempty(oletter)                                    
               dwelling = dwelling+1;                    
               oletter = trl.letter;
           elseif strcmp(oletter,trl.letter) 
               dwelling = dwelling+1;                    
               oletter = trl.letter;
           else
               oletter = trl.letter;
               dwelling = 0;
           end   

           if dwelling>=trl.DwellingCriterion

                trl.word = char(strcat({trl.word},{trl.letter}));  
                dwelling = 0;
                EVT = trl.letter;

           end
       end       
       POS.EVT(1:end-1) = POS.EVT(2:end);
       POS.EVT(end) = EVT;
 
    end
    function nDrawFeedback(win,gfx,trl)

        csz = dwelling + 1;
        minsize = 5;
        Screen('DrawDots', win, [trl.fx, trl.fy], fix(10.*(csz/trl.DwellingCriterion))+minsize, trl.ccolor);
        Screen('TextSize', win, 15);
        if ~isempty(layerid) && ~isempty(letterid) && (dwelling > trl.DwellingCriterion*(2/3)) % show the letter when close to criterion
            DrawFormattedText(win, gfx.Let{layerid}(letterid), 'center', 'center', ...
                           [255 255 255],[],[],[],[],[],...
                           CenterRectOnPoint([0 0 10 10], trl.fx, trl.fy)); 
        end     
        if ~isempty(trl.word)        
            Screen('TextSize', win, 30);
            DrawFormattedText(win, trl.word, 'center', gfx.wrdCentre(2), [255 255 255], 80, [],[],2);            
        end
        Screen('DrawingFinished', win);
        vbl = Screen('Flip', win, vbl + 0.5 * gfx.ifi, 1);        
    end
    function trl=nDisplayWord(win,gfx,trl)
        trl.word_to_write = trl.sWORDS{trl.trl};        
        Screen('TextSize', win, 30);
        DrawFormattedText(win, trl.word_to_write, 'center', gfx.wrdCentre(2), [255 255 255],80,[],[],2);
        vbl = Screen('Flip', win, vbl + 0.5 * gfx.ifi);        
        WaitSecs(.5);                                                      % minimum of half a second, press key to go on 
    end
    function trl = nGetResp(trl,el,util)
        [ ~, ~, keyCode ] = KbCheck;
        if numel(trl.word)==6
            if strcmp(trl.word,trl.word_to_write)
                trl.DwellingCriterion=trl.DwellingCriterion-2;
                fprintf('DwellingCrit %d\n', trl.DwellingCriterion);
                Beeper(800,.5,0.1) % correct, jaunty tone :-)
            else
                trl.DwellingCriterion=trl.DwellingCriterion+2;
                fprintf('DwellingCrit %d\n', trl.DwellingCriterion);
                Beeper(200,1,0.1)  % wrong, sad tone :-(
            end
            trl.complete = 1;
        elseif keyCode(KbName('c'))
            if(trl.Eyelink), EyelinkDoTrackerSetup(el); end         % Calibrate and Validate               
        elseif keyCode(KbName('q'))
            fCloseEL(util);            
            sca;
            clear;
        end
    end

    % eye movements / filter
    function trl = nDetectFix(trl,el)    
        if Eyelink('NewFloatSampleAvailable') > 0            
            evt = Eyelink('newestfloatsample');

%        time: 4173808
%        type: 200
%       flags: 26369
%          px: [-32768 -32768]
%          py: [-32768 -32768]
%          hx: [-32768 -32768]
%          hy: [-32768 -32768]
%          pa: [32768 0]
%          gx: [-32768 -32768]
%          gy: [-32768 -32768]
%          rx: 52.6000
%          ry: 70.1000
%      status: 32768
%       input: 32768
%     buttons: 0
%       htype: 0
%       hdata: [0 0 0 0 0 0 0 0]
            
            x = evt.gx(el.eye_used+1);
            y = evt.gy(el.eye_used+1);     

            if x~=el.MISSING_DATA && y~=el.MISSING_DATA    
                trl.x = x;
                trl.y = y;                         
            else
                trl.x = ox;
                trl.y = oy;
            end            

        end
        %fprintf('trl.x %4.2f, trl.y %4.2f\n', trl.x, trl.y);
    end
end