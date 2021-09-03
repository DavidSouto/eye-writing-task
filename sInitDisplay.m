%% INIT DISPLAY & GAMMA CORRECTION
if trl.DEBUG
    Screen('Preference','SkipSyncTests', 1);
    res=[0 0 1280 1024];
else
    Screen('Preference','SkipSyncTests', 0);
end

AssertOpenGL;

PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
        
screenID = max(Screen('Screens'));

%keyboard
%Screen('OpenWindow',windowPtrOrScreenNumber [,color] [,rect][,pixelSize][,numberOfBuffers][,stereomode][,multisample][,imagingmode][,specialFlags][,clientRect]);
[win winRect] = Screen('OpenWindow', screenID, [128 128 128],res,32,2,[],[],[],[],[]);

Screen(win,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

HideCursor;    
KbName('UnifyKeyNames');

Screen('TextFont', win, 'Arial');
Screen('TextSize', win, 20);