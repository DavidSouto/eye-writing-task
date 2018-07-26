function [newGamma oldGamma] = exp_mon_init(win,epar)
    %initMon(epar.GAMMA_TABLE);
    if epar.doGamma,
    	newGamma(:,1) = dlmread([epar.GAMMA_TABLE '.r']);
    	newGamma(:,2) = dlmread([epar.GAMMA_TABLE '.g']);
    	newGamma(:,3) = dlmread([epar.GAMMA_TABLE '.b']);
    	newGamma = newGamma./255;
    	oldGamma = Screen('LoadNormalizedGammaTable',win,newGamma);
    end
end