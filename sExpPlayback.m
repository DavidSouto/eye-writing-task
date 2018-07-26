%% (EL2_playback: courtesy of Alexander C Schutz @ Giessen) %%%%%%%%%%%
if ~trl.noEyelink
    %['\data\e',int2str(trl.exp),'s',int2str(trl.sub),
    %util.el_datname]
    status = EL2_playback(util.el_datdir);
    if ~status
        error('playback failed');
    else
        util.el_datdir
        eyed = textread(util.el_datdir);

        errcodes.SAC = 1;
        errcodes.BLINK = 2;
        errcodes.MSG = 3;
        errcodes.SAC_AMP = 4;
        errcodes.GAIN = 5;

        % SAC index
        Z = bitget(eyed(:,4),errcodes.SAC);
        SAC = bitget(eyed(:,4),errcodes.SAC);
        BLINK = bitget(eyed(:,4),errcodes.BLINK);                
        EVE = find(bitget(eyed(:,4),errcodes.MSG)==1);

        SACCADE = [];
        if ~isempty(SAC), 
            % old (does not work (?))
            %SACCADE(:,1) = SAC(abs(diff([-1000;SAC]))>1);
            %SACCADE(:,2) = SAC(abs(diff([SAC;10000]))>1);

            % beginning (wrong onset at start)
            SACCADE(:,1) = find(diff([0;SAC]) == 1);

            % end (wrong offset at end)
            SACCADE(:,2) = find(diff([SAC;0]) == -1);   

        end

        % gaze (x,y)
        X = eyed(:,2);
        Y = eyed(:,3);

        % center and convert to degrees
        X = (X-gfx.scMid(1))./gfx.h_deg2pix;
        Y = -(Y-gfx.scMid(2))./gfx.v_deg2pix;     

        % get saccade amplitudes (2D amp)
        if ~isempty(SACCADE),
            % GET X AMPLITUDE (X1-X0)
            SACCADE(:,3) = X(SACCADE(:,2))-X(SACCADE(:,1));
            % GET Y AMPLITUDE (Y1-X0)
            SACCADE(:,4) = Y(SACCADE(:,2))-Y(SACCADE(:,1));
            % NORM: 
            SACCADE(:,5) = sqrt(SACCADE(:,3).^2 + SACCADE(:,4).^2);
        end

        % (assuming 1000 hz)
        exp.saccade_crit = 1;

        exp.start = EVE(2)-200/(1000/gfx.sfreq); % start quality check from cohere onset and a little before
        exp.end = EVE(3);

        % CHECK SACC, BLINKS AND AMP DURING CRITICAL INTERVAL
        if sum(BLINK(exp.start:exp.end)), 
            fprintf('blink error during critical interval\n');
            trl.err = bitset(trl.err,errcodes.BLINK); 
        end
        if sum((SAC > exp.start) & (SAC < exp.end)), % if there is a saccade within (or straddling) the interval
            for n=1:length(SACCADE(:,1)),            % saccades length                        
                fprintf('%.2f degrees sacc (crit %.2f)\n',SACCADE(n,5), exp.saccade_crit);
                if SACCADE(n,5)>exp.saccade_crit,
                    trl.err = bitset(trl.err,errcode.SAC_AMP); 
                end
            end
        end

        %------- Do the on-line pursuit quality check
        % add a security margin (20 samples) to non-saccade
        % bits = NaN
        DX = diff(X).*gfx.sfreq;
        [b,a] = butter(2,40/(gfx.sfreq/2),'low');
        DX  = filtfilt(b, a, DX);
        if sum(Z), % if at least one saccade episode
            SSAC = sexpandsacc(Z,SACCADE,20);
            DX(logical(SSAC)) = NaN; % no interpolation, traces marked as NaN for the extended duration of the saccade
        end
        Gain = nanmean(DX(exp.start:exp.end)./gfx.vel);
        if ~trl.fixation,
            if abs(Gain) < 0.8,
                trl.err = bitset(trl.err, errcodes.GAIN);
                fprintf('gain too low: %.2f\n', abs(Gain));
            end
            trl.rs = abs(nanmean(DX(exp.start:exp.end))) - gfx.vel;
            trl.rs
            trl.gain = Gain;
            nanmean(DX(exp.start:exp.end))
        end
    end                

end