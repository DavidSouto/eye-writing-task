function [trl, POS]= fDrawStimProc(win,gfx,trl,el)
   
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
    POS.evt = NaN*ones(5000,1); % get character with char()
    
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
        Eyelink('Message','SYNCTIME'); % EVENT 2
    end
    
    priorityLevel=MaxPriority(win);
    Priority(priorityLevel);

	trl.complete = 0;
    trl.repeat = 0;
    trl.word = '';
    SetMouse(gfx.wrdCentre(1),gfx.wrdCentre(2));                           % Init mouse location at center
    
    try
        nDisplayKeyboard(win,gfx);                                         % dial virtual keyboard
    catch me
        sca
        keyboard
    end
    trl = nDisplayWord(win,gfx,trl);     % display word to write
    
    if trl.joystickControl
        h(5) = 0;
        while(~h(5))
            h = jst;        
        end
    else
        KbWait([], 2);
    end
    
    % TBD: Need drift correction routine that does not cause a redraw
    %result=EyelinkDoDriftCorrect(el, gfx.scMid(1),gfx.scMid(2))
    
    evt = '!'; % trial start event
    t0 = GetSecs();
    [trl POS]=nGetCursorCoords(trl,evt,POS);                               % cursor x,y coordinates, from eye or mouse    
    [trl evt]=nTestSelection(gfx,trl,evt);    
    while (GetSecs()-t0)<0.5 % continue updating cursor position 1 second after complete
        evt = '0';
        [trl evt]=nTestSelection(gfx,trl,evt);
        [trl POS]=nGetCursorCoords(trl,evt,POS);                           % cursor x,y coordinates, from eye or mouse
        nDrawFeedback(win,gfx,trl); % show cursor position
        if ~trl.complete
            t0 = GetSecs();
            trl=nGetResp(trl,el);                                             % opportunity to skip if incomplete
        end  
    end
    fprintf('trial %d repeat %d\n', trl.trl, trl.repeat);
    if trl.Eyelink
        Eyelink('Message','TRIAL_END');  % EVENT 3
        Eyelink('command','end_realtime_mode()'); 
    end
    Priority(0);
    
    % nested functions 
    function [trl, POS] = nGetCursorCoords(trl,evt,POS)
        if trl.control                                                     % eye-control mode
            trl = nDetectFix(trl,el);  
        else                                                               % mouse-control mode
            [x, y, ~] = GetMouse; 
            trl.x = x;            
            trl.y = y;                        
        end
                
        %% --- Exponential filtering ---
        % TBD: DX(logical(sacc_id)) = NaN;
        % TBD: start with a simulation 
        % alpha = 0.4;
        
        %lseries = length(DX);
        %l = ones(lseries,1);
        %err = ones(lseries,1);
        
        %l(1) = DX(3)-DX(2);
        %for t = 3:length(DX)
        %    if isnan(DX(t))
        %       l(t) = l(t-1); 
        %    else
        %       err(t) = DX(t) - l(t-1); 
        %
        %       l(t) = l(t-1) + alpha * err(t);
        %    end
        %end                    
        
        %% --- Update position buffer with last sample ---
        try
            POS.x(1:end-1) = POS.x(2:end);                                     
            POS.y(1:end-1) = POS.y(2:end);   
            POS.t(1:end-1) = POS.t(2:end);        
            POS.evt(1:end-1) = POS.evt(2:end);        
            POS.x(end) = trl.x;                                                % A 5000 element BUFFER RECORDING CURSOR MOVEMENTS
            POS.y(end) = trl.y;
            POS.t(end) = GetSecs();
            POS.evt(end) = evt;
        catch
            sca
            keyboard
        end
    end
    function nDisplayKeyboard(win,gfx)
        LPix = 30;
        Screen('TextSize', win, LPix);
 
        for ly = 1:4                                                                 % layer

            Screen('FrameArc',win, 0, CenterRectOnPoint([0 0 gfx.bE(ly+1)*2 gfx.bE(ly+1)*2], gfx.kbCentre(1), gfx.kbCentre(2)),...
                    (gfx.b(ly).C(1)-gfx.bDa(ly)*.5).*180/pi, ...% start
                    (gfx.maxAngle(ly)+gfx.bDa(ly))*180/pi,3);                        % arc

            Screen('FrameArc',win, 0, CenterRectOnPoint([0 0 gfx.bE(ly)*2 gfx.bE(ly)*2], gfx.kbCentre(1), gfx.kbCentre(2)),...
                    (gfx.b(ly).C(1)-gfx.bDa(ly)*.5).*180/pi, ...% start
                    (gfx.maxAngle(ly)+gfx.bDa(ly)).*180/pi,3);                       % arc

            % need to rotate relative to FrameArc coords
            [fH, fV]=pol2cart([gfx.b(ly).C-gfx.bDa(ly)*.5-pi/2, gfx.b(ly).C(end)+gfx.bDa(ly)*.5-pi/2], gfx.bE(ly) ); % TBD: one is missing ...
            [tH, tV]=pol2cart([gfx.b(ly).C-gfx.bDa(ly)*.5-pi/2, gfx.b(ly).C(end)+gfx.bDa(ly)*.5-pi/2], gfx.bE(ly+1));

            [letterx,lettery]=pol2cart(gfx.b(ly).C-pi/2, gfx.bE(ly)+gfx.bD(ly)*.5);                                                            

            for n = 1:length(gfx.Let{ly})                                               % Display letters
               Screen('DrawLines', win, [[fH(n);fV(n)], [tH(n);tV(n)]], 3, 0, gfx.kbCentre,2);        

               DrawFormattedText(win, gfx.Let{ly}(n), 'center', 'center', ...
                   [255 255 255],[],[],[],[],[],...
                   CenterRectOnPoint([0 0 LPix LPix], letterx(n)+gfx.kbCentre(1), lettery(n)+gfx.kbCentre(2)));   
            end            
            Screen('DrawLines', win, [[fH(end);fV(end)], [tH(end);tV(end)]], 3, 0, gfx.kbCentre,2);        

        end
        vbl = Screen('Flip', win, vbl + 0.5 * gfx.ifi, 1);                 
        %KbWait([],2);

    end
    function [trl evt] = nTestSelection(gfx,trl,evt)
        layerid =[];        
        letterid=[];
        [th, r] = cart2pol(trl.x-gfx.kbCentre(1),trl.y-gfx.kbCentre(2));   % Polar coordinates
        if (r<gfx.bE(end)) && (r>gfx.bE(1))                                % Detect layer
                                                                
            layerid=sum(r>gfx.bE);
            
            % rotate
            minLetterAngle=utilPangle(gfx.b(layerid).C(1)-gfx.bDa(layerid)*.5-pi/2);
            maxLetterAngle=utilPangle(gfx.b(layerid).C(end)+gfx.bDa(layerid)*.5-pi/2);
                            
            pangle = utilPangle(th);
            %% All positive angles
            if ~(pangle< minLetterAngle) && ~(pangle > maxLetterAngle)     % Detect letter

                   cangle = utilPangle(gfx.b(layerid).C-pi/2);

                   AngleDiff = abs(cangle-pangle);                                  
                   minAngleDiff = AngleDiff==min(AngleDiff);               % find the letter closest to centroid 
                   letterid = minAngleDiff;   
                   trl.letter = gfx.Let{layerid}(letterid);

                   if  isempty(oletter)                                    % Dwell time criterion
                       dwelling = dwelling+1;                    
                       oletter = trl.letter;
                   elseif strcmp(oletter,trl.letter) 
                       dwelling = dwelling+1;                    
                       oletter = trl.letter;
                   else
                       oletter = trl.letter;
                       dwelling = 0;
                   end   

                   if strcmp(trl.letter,'<') && dwelling>=trl.DwellingCriterion       % Letter selection

                        if length(trl.word)<=1
                            trl.word = [];
                        else
                            trl.word=trl.word(1:end-1);
                        end
                        dwelling = 0;
                        evt = trl.letter;

                   elseif strcmp(trl.letter,'_') && dwelling>=trl.DwellingCriterion

                        trl.word = char(strcat({trl.word},{' '}));  
                        dwelling = 0; 
                        evt = trl.letter;

                   elseif dwelling>=trl.DwellingCriterion
                       
                        trl.word = char(strcat({trl.word},{trl.letter}));  
                        dwelling = 0;
                        evt = trl.letter;

                   end                 
            end            
        end
    end
    function nDrawFeedback(win,gfx,trl)

        csz = dwelling + 1;
        minsize = 5;
        Screen('DrawDots', win, [trl.x, trl.y], fix(10.*(csz/trl.DwellingCriterion))+minsize, trl.ccolor);
        Screen('TextSize', win, 15);
        if ~isempty(layerid) && ~isempty(letterid) && (dwelling > trl.DwellingCriterion*(2/3)) % show the letter when close to criterion
            DrawFormattedText(win, gfx.Let{layerid}(letterid), 'center', 'center', ...
                           [255 255 255],[],[],[],[],[],...
                           CenterRectOnPoint([0 0 10 10], trl.x, trl.y)); 
        end                   
        if ~isempty(trl.word)        
            Screen('TextSize', win, 30);
            DrawFormattedText(win, trl.word, 'center', gfx.wrdCentre(2), [255 255 255], 80, [],[],2);            
        end
        Screen('DrawingFinished', win);
        vbl = Screen('Flip', win, vbl + 0.5 * gfx.ifi, 1);        
    end
    function trl=nDisplayWord(win,gfx,trl)        
        trl.word_to_write = trl.sWORDS{trl.session_corpus_id(trl.trl)};        
        %Screen('TextColor',win,[255 0 0]);
        Screen('TextSize', win, 30);
        DrawFormattedText(win, trl.word_to_write, 'center', gfx.wrdCentre(2), [255 255 255],80,[],[],2);
        vbl = Screen('Flip', win, vbl + 0.5 * gfx.ifi);        
        WaitSecs(.5);                                                      % minimum of half a second, press key to go on 
    end
    function trl = nGetResp(trl,el)
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
        end            
%         if keyCode(KbName('RightArrow'))
%             trl.complete = 1; % go to next word
%         elseif keyCode(KbName('LeftArrow'))
%             trl.repeat = 1; % go back: show the same sequence agai
%         elseif keyCode(KbName('UpArrow'))
%             trl.DwellingCriterion = trl.DwellingCriterion + 5;
%             Screen('TextSize', win, 30);
%             dtext = 'increased dwelling time';
%             DrawFormattedText(win, dtext, 'center', gfx.wrdCentre(2)-50, [255 255 255], 80, [],[],2);            
%         elseif keyCode(KbName('DownArrow'))
%             trl.DwellingCriterion = trl.DwellingCriterion - 5;
%             if trl.DwellingCriterion <=0
%                 trl.DwellingCriterion = 5;
%             end
%            Screen('TextSize', win, 30);
%            dtext = 'decreaset dwelling time';
%            DrawFormattedText(win, dtext, 'center', gfx.wrdCentre(2)-50, [255 255 255], 80, [],[],2);            
        KbName('UnifyKeyNames');
        [~, ~, keyCode, ~] = KbCheck;
        if keyCode(KbName('q'))
            sca
            clear;
        elseif keyCode(KbName('c'));
            
            if(trl.Eyelink), EyelinkDoTrackerSetup(el); end         % Calibrate and Validate   
 
        end
    end
    % eye movements / filter
    function trl = nDetectFix(trl,el)    
        if Eyelink('NewFloatSampleAvailable') > 0
            
            evt = Eyelink('newestfloatsample');
%sca
%keyboard
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
            %Col.xscaled = 15; % scaled signal
            %Col.yscaled = 17; % scaled signal
            
            if x~=el.MISSING_DATA && y~=el.MISSING_DATA
                                
               % we want position in the screen coordinates for once
                trl.x = x;%./gfx.h_deg2pix;
                trl.y = y;%-gfx.scMid(2));%./gfx.v_deg2pix;                
            
                %trl.t_new = evt.time;                                
            else
                %trl.x = NaN; trl.y = NaN; trl.t_new = NaN;                
                % keep the old position
                trl.x = ox;
                trl.y = oy;
%                 if isnan(ox);
%                     trl.x = x;
%                     trl.y = y;
%                 end
            end            
%             l(1) = DX(3)-DX(2);
%             for t = 3:length(DX)
%                 if isnan(DX(t))
%                    l(t) = l(t-1); 
%                 else
%                    err(t) = DX(t) - l(t-1); 
% 
%                    l(t) = l(t-1) + alpha * err(t);
%                 end
%             end
%             if isnan(ox);
%                 ox = trl.x;
%                 oy = trl.y;
%             end
%             trl.x = trl.alpha * trl.x + (1-trl.alpha)*ox;
%             trl.y = trl.alpha * trl.y + (1-trl.alpha)*oy;
%trl.ox = trl.x;
%            trl.oy = trl.y;

        else
%            trl.x = trl.ox;
%            trl.y = trl.oy;
        end
        fprintf('trl.x %4.2f, trl.y %4.2f\n', trl.x, trl.y);

    end
end