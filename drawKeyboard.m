% TBD: DRAW THE KEYBOARD 
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