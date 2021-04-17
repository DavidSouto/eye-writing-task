% Feb 2018 DS @ Uni Leicester
%
% USAGE:
% controlMode = 0 (mouse), 1 (eye), Eyelink = 0 (no eyelink), 1(eyelink)

% DESCRIPTION: 
% - calibration switch during experiment (c)

% gamma correction is off
function fExp(subject, session, control_mode)

    if(nargin < 3), error('missing arguments'); end
    
    trl.exp = 1; 
    trl.sub = subject;
    trl.blk = session;
    trl.control = control_mode;    
    trl.Eyelink = 1; 
    trl.DEBUG = 1;
    trl.joystickControl = 0;

    trl.DwellingCriterion = 30;                             % Dwelling criterion 
    
    if trl.control== 1 && ~trl.Eyelink
        error('Wrong parameter, Eyelink has to be on to do the eye-writing');
    end

    if trl.joystickControl
        trl.joy = vrjoystick(1);
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
            while( k < trl.nTrials)
                trl.trl = k;
                
                if(trl.Eyelink), sStartRec; end                            % start eyelink rec
                                                                           % display stimuli
                %EyelinkDoDriftCorrection(el); %see how to control
                [trl, POS]= fDrawStimProc(win,gfx,trl,el);

                if(trl.Eyelink), sStopRec; end                             % stop eyelink rec                             
                if(trl.Eyelink), forceCalib(el,trl); end                   % allow to force calibration if necessary
                
                trl = fSaveTrial(trl, util, POS);                          
                if ~trl.repeat
                    k = k + 1;                                             % only move on in that case, but store an eye link file: trl.num
                end
            end  
    end
    if trl.Eyelink, fCloseEL(util); end    
    fQuitProg();
end
function trl=fLoadCorpus(trl)                                              % create corpus
    fid=fopen('ucrel.corpus.csv');
    WORDS = textscan(fid,'%s %d %s','delimiter',',');
    fclose(fid);

    % make a random selection of trl.ntrl words over the whole n-letter
    % corpus, repetition is allowed, should make an equal distribution of
    % frequencies accross sessions (combination with replacement)
    id = find(WORDS{2}==6); % 765 words
           
    for n = 20:(200+2*trl.TrainingTrials+20) %length(id)                                                   
        trl.sWORDS{n-19}= WORDS{3}{id(n)};% word shortlist
    end
    % for that block
    if trl.training
        trl.session_corpus_id = 1:trl.TrainingTrials;
    else       
        if trl.blk>0
            nsession = 1;
        else
            nsession = 0;
        end
        trl.session_corpus_id = (trl.TrainingTrials*nsession+(nsession-1)*trl.nTrials+1):(trl.TrainingTrials*nsession+nsession*trl.nTrials); 
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
    
    util.clob_name = sprintf('%sc.%ds%d.mat',util.dir_name, trl.blk, trl.num);

    % save cursor traj
    DATA = [POS.x, POS.y, POS.t, POS.evt];
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
                    '---------------\n Press key to proceed\n'];

    DrawFormattedText(win, screen_param, gfx.scMid(1)-gfx.scMid(1)*.5, gfx.scMid(2)-400, 255,[],[],[],2); 
    Screen('Flip',  win, [], 1);
    
    WaitSecs(.5);
    
    if trl.joystickControl
        
        %h(5) = 0;
        h = 0;
        while(~h)%(h(5)))
            h = button(trl.joy,1);
            %h = jst;        
        end
    else
        KbWait([], 2);
    end
end
function fDispBlockInfo(win,trl,expon)
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

    WaitSecs(1);
    
    if trl.joystickControl   
        %h(5) = 0;
        h = 0;
        while(~h)%(h(5)))
            h = button(trl.joy,1);
            %h = jst;        
        end
    else
        KbWait([], 2);
    end

end
function fCloseEL(util)
    Eyelink('command', 'generate_default_targets = YES');
    Eyelink('CloseFile');
    Eyelink('ReceiveFile',[],util.dir_name,1);  
    Eyelink('Shutdown');
    sca;
    commandwindow; 
end
function forceCalib(el,trl)
        [~,~,KbCode] = KbCheck;
        if KbCode(KbName('b')) 
            if(trl.Eyelink)
                EyelinkDoTrackerSetup(el);                                 % Calibrate & Validate                
            end
        elseif KbCode(KbName('q'))                                         % Allow early exit
            sca;
            clear;
        end
end 
