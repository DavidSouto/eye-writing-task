function [trl] = fGetResp(trl)

        [dummy KbCode] = KbWait([], 2);  

        if KbCode(KbName('LeftArrow')) 

            trl.resp = 0;
            
        elseif KbCode(KbName('RightArrow')),

            trl.resp = 1;
        
        elseif KbCode(KbName('q'));     
            
            trl.resp = NaN;
            sca;
            clear all;

        elseif KbCode(KbName('space')),
           
            trl.loop = 0;
            
            return;
            
        else
            
            trl.resp = NaN;
            
        end
end    