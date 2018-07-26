% CUSTOMCALIBRATION for eye tracking using EyeLink
%
% This code should show you how to get a custom calibration in real time.
% Stick this piece of code in your eyetracking setup. The way it works is
% you pick the number of dots (samples), and then move the mouse around
% (there will be no cursor, just a simulation of what the calibration dot
% will look like) until you find a good place for the dot during the
% calibration/validation.
% 
% A good location will be determined by where the subject can actually see
% on the screen (in case part of it is occluded) and general knowledge
% about what makes a good calibration (do the corners and one in the
% middle). When you found a good spot, click the mouse.
% 
% There will be a counter in the top left telling you how many to go.
% 
% ADVICE: Because the EyelinkDoTrackerSetup comes after all of this, make
% sure that you have the eye tracker all set up so that it is physically
% oriented correctly and in focus so that the calibration dots are in the
% right location for the tracker.
%
%
% EG Gaffin-Cahn
% 06/2014

% ... other eyelink commands and updating of defaults

%Eyelink('Command', 'generate_default_targets = NO');

% instructions: move the mouse around and talk to subject until you find a
% good place. When you have a good place for a calibration spot, click the
% mouse.
nsamples = 13;
sample = 1;
B = 20;
S = 5;
drawKeyboard;
xc = NaN.*ones(1,nsamples);
yc = NaN.*ones(1,nsamples);
[x,y,buttons] = GetMouse(win);
drawKeyboard;
while sample <=nsamples; 
        
    [x,y,buttons] = GetMouse(win); 

    drawKeyboard;
    Screen('DrawText',win,sprintf('%d / %d',sample,nsamples),50,50);
    Screen('FillOval',win,[0 0 0],[x-B, y-B, x+B, y+B]);
    Screen('FillOval',win,[128 128 128],[x-S, y-S, x+S, y+S]);
    Screen('Flip',win);

    if any(buttons)
        xc(sample) = x; yc(sample) = y;
        sample = sample + 1;
        WaitSecs(0.5);
    end

end
%sca
%keyboard
Eyelink('Command','calibration_samples = %d',nsamples);
Eyelink('Command','calibration_sequence = 0,1,2,3,4,5,6,7,8,9,10,11,12');
Eyelink('Command','calibration_targets = %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',...
    xc(1),yc(1),  xc(2),yc(2),  xc(3),yc(3),  xc(4),yc(4),  xc(5),yc(5), xc(6),yc(6),  xc(7),yc(7), xc(8),yc(8),  xc(9),yc(9),xc(10),yc(10), xc(11),yc(11),xc(12),yc(12),xc(13),yc(13));
Eyelink('Command','validation_samples = %d',nsamples);
Eyelink('Command','validation_sequence = 0,1,2,3,4,5,6,7,8,9,10,11,12');
Eyelink('Command','validation_targets = %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',...
    xc(1),yc(1),  xc(2),yc(2),  xc(3),yc(3),  xc(4),yc(4),  xc(5),yc(5), xc(6),yc(6),  xc(7),yc(7), xc(8),yc(8),  xc(9),yc(9),xc(10),yc(10), xc(11),yc(11),xc(12),yc(12),xc(13),yc(13));

% other eyelink commands and updating of defaults ...