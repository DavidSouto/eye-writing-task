% August 2018 DS @ Uni Leicester
%
% USAGE:
% controlMode = 0 (mouse), 1 (eye), Eyelink = 0 (no eyelink), 1(eyelink)

% DESCRIPTION: 
% - calibration switch during experiment (c)
% - Instructions: follow the shortest path, fill skipped letters with _
% Optional 13 point calibration, if calib = 1

% gamma correction is off
function fExp(subject, session, control_mode)

    trl.DEBUG = 1; % ZERO FOR EXP
    trl.joystickControl = 0;
    trl.Eyelink = 0; 
    
    if(nargin < 3), error('missing arguments'); end
    
    trl.exp = 1; 
    trl.sub = subject;                                                                                                                                                      
    trl.blk = session;
    trl.control = control_mode;    
    if ~trl.Eyelink
        warning('DEMO: NO EYELINK MODE');
    end
    
    %% ADJUST CALIB SIZE HERE
    gfx.CALIB_X = 10;
    gfx.CALIB_Y = 5;    
    
    if nargin<4
                                                                                                                                            
        trl.calib = 1;
    else
        trl.calib = calib; 
    end

    trl.criterion_filename = [int2str(trl.sub) '.savecrit.mat'];% Dwelling criterion 
    if exist(trl.criterion_filename,'file') % load the last value of previous session if exists
        load(trl.criterion_filename,'VAL');
        trl.DwellingCriterion = VAL;
    else
        trl.DwellingCriterion = 40;                         % Dwelling criterion 
    end

    sInitDisplay;                                           % Initialize display
    sExpVars;                                               % Experimental (trl) and Display variables (gfx)
        
    if(trl.Eyelink), sInitEL; else el = []; end             % Initialize eyelink 
    if(trl.Eyelink), EyelinkDoTrackerSetup(el); end         % Calibrate and Validate   
        
    fDispExpInfo(win,trl,gfx);    
    for expon = 0:1
            trl.training = 1-expon;
            if(expon==0) 
                trl.nTrials = trl.TrainingTrials; 
            elseif(expon==1) 
                trl.nTrials=trl.ntrl; 
            end
            trl=fLoadCorpus(trl);
            trl.expon = expon;
            fDispBlockInfo(win, trl, expon);
            k = 1;
            while( k <= trl.nTrials)
                trl.trl = k;
                if(trl.Eyelink)
                    if trl.joystickControl
                        el.allowjoystick = 1;
                    else
                        el.allowjoystick = 0;
                    end
                    EyelinkDoDriftCorrectionDS(el,gfx.wrdCentre(1),gfx.scMid(2)); %see how to control
                end
                if(trl.Eyelink), sStartRec; end                            % start eyelink rec
                                                                           % display stimuli
                [trl, POS]= fDrawStimProc(win,gfx,trl,util,el);

                if(trl.Eyelink), sStopRec; end                             % stop eyelink rec                             
                if(trl.Eyelink), forceCalib(el,trl,util); end                   % allow to force calibration if necessary
                
                trl = fSaveTrial(trl, util, POS);                          
                if ~trl.repeat
                    k = k + 1;                                             % only move on in that case, but store an eye link file: trl.num
                end
            end  
    end
    if trl.Eyelink, fCloseEL(util); end    
    fQuitProg();
end
function trl=fLoadCorpus(trl)     % create the corpus

    if IsLinux
        fid=fopen('corpus/ucrel.corpus.filt.csv');
    else
        fid=fopen('corpus\ucrel.corpus.filt.csv');
    end
    WORDS = textscan(fid,'%s %d %d %2.1f %d','delimiter',',');
    fclose(fid);

    % make a random selection of trl.ntrl words over the whole n-letter
    % corpus, repetition is allowed, should make an equal distribution of
    % frequencies accross sessions (combination with replacement)
    id = find(WORDS{5}==6); % 765 6-letter words
    word_sample = randperm(length(id),trl.NWORDS);
    for n = 1:length(word_sample)
        trl.sWORDS{n} = WORDS{1}{id(word_sample(n))}; % sampling without replacement
    end
end
function fQuitProg()
    sca;
    clear;
end
function [trl] = fSaveTrial(trl, util, POS)

    fid = fopen(util.log_name,'a');
    fprintf(fid, '%d %d %d %d %d %d %s %s %d %d %d\n', ... 
                trl.exp, trl.expon, trl.sub, trl.blk, trl.num, trl.trl, trl.word, trl.word_to_write, trl.repeat, trl.DwellingCriterion, trl.control);     
    fclose(fid);
    
    % trl.session_corpus_id = (trl.TrainingTrials*nsession+(nsession-1)*trl.nTrials+1):(trl.TrainingTrials*nsession+nsession*trl.nTrials);    
    util.clob_name = sprintf('%s/c.%ds%d.mat',util.dir_name, trl.blk, trl.num);

    DATA = [POS.x, POS.y, POS.t, POS.EVT];
    save(util.clob_name,'DATA');     
            
    trl.num = trl.num + 1;
end
function fDispExpInfo(win,trl,gfx)   

    Screen('Flip',  win, [], 1);
    Screen('TextFont', win, 'Arial');
    Screen('TextSize', win, 20);
    screen_param = ['---------------\n'...
                    'Exp: ' int2str(trl.exp) '\n'...
                    'Monitor frequency: ' int2str(gfx.monitor_freq) ' Hertz\n'...
                    'Screen resolution: ' int2str(gfx.h_pixel) 'x' int2str(gfx.v_pixel) ' pixels\n'...
                    'Desired sub. distance: ' int2str(gfx.sub_distance) ' cm\n'...
                    'Designed for screen: ' gfx.screen '\n' ...
                    '---------------\n Press A on the gamepad to continue\n'];

    DrawFormattedText(win, screen_param, gfx.scMid(1)-gfx.scMid(1)*.5, gfx.scMid(2)-400, 255,[],[],[],2); 
    Screen('Flip',  win, [], 1);
    
    WaitSecs(.5);
    
    if trl.joystickControl
        fWaitForJoyStick()
    else
        KbWait([], 2);
    end
end
function fDispBlockInfo(win,trl,expon)
    
    if trl.joystickControl   
        fWaitForJoyStick()
    else
        KbWait([], 2);
    end
    Screen('Flip',  win); % clear the screen
    Screen('TextSize', win, 30);

    PHASE = {'TRAINING PHASE: ', 'EXPERIMENTAL PHASE: '};

    if trl.control
        INSTRUCTION = 'your task is to write target\n by looking at letter boxes\n \n Selection is made by dwelling time\n Dwelling time is reduced when you are successful\n and increased when you are unsuccessful\n';
    else
        INSTRUCTION = 'your task is to write target words\n by moving the mouse over letter boxes\n Selection is made by dwelling time\n Dwelling time is reduced when you are successful\n and increased when you are unsuccessful\n'; 
    end
    DrawFormattedText(win, sprintf('%s%s\n', PHASE{expon+1}, INSTRUCTION), 'center', 'center', [255 255 255],80,[],[],2);
    Screen('Flip',  win, [], 1);

    WaitSecs(5);

    if trl.joystickControl   
        fWaitForJoyStick()
    else
        KbWait([], 2);
    end
end
function forceCalib(el,trl,util)
        [~,~,KbCode] = KbCheck;
        if KbCode(KbName('b')) 
            if(trl.Eyelink)
                EyelinkDoTrackerSetup(el);                                 % Calibrate & Validate                
            end
        elseif KbCode(KbName('q'))                                         % Allow early exit
            fCloseEL(util);
            sca;
            clear;
        end
end
