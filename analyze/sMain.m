%% Init 
sDefine;

Trl.match = 1;

lab = 0;
if lab
    addpath(genpath('/media/dsadmin/My Files1/src/MATLAB/'));
else
    addpath(genpath('c:/Users/Asus/Documents/src/MATLAB/'));
end

for subcnt = 1:length(sub_list) 

    S = struct([]);

    Trl.sub    = sub_list(subcnt); % subject number
    datafile = [int2str( sub_list(subcnt) ) '.mat'];

    if ~exist(datafile,'file') % manually delete the file to re-analyze

        Man.log = sprintf('%se%ds%d/e%ds%d.log', Man.datpath, Trl.exp, Trl.sub, Trl.exp, Trl.sub);

        DATA = textread(Man.log);

        MDATA = DATA(DATA(:,Col.train)==0 & DATA(:,Col.task)==0,:); 
        DATA = DATA(DATA(:,Col.train)==0 & DATA(:,Col.task)==1,:);        

        [nrows, ~] = size(DATA);

        % init trace struct
        for cc = 1:3
            for tt = 1:3
                S(cc,tt).cnt = 0;
            end
        end
        
        % FIND PERCEPTUAL MATCHES
        sPercMatches;
        for counter=1:nrows

            Trl.num = counter;
            Trl.trl = Trl.num;
            error(Trl.num) = 0;       

            sLoad; 

            sAnalyze;    % do pursuit/fixation quality parsing
            if (max(abs(intdx))>100 || max(abs(intdy))>100) && error(Trl.num)==0
                
                figure; plot(intdx);
                
            end
            
            S(cid,speedid).cnt = cond_inc + 1;
    
            timetrl = stimOff-stimOn; %length(X); % detetect wrong timestamp
            fprintf('trial %d processed, condition %d, subject %d error: %d length trial %d\n',Trl.num, DATA(Trl.num,Col.cond), Trl.sub, error(Trl.num), timetrl);

            for n = 1:length(Err.TYPES) % go through errors
                if bitget(error(Trl.num),n)
                    fprintf('error %s\n', char(Err.TYPES(n)));
                end
            end
        end
        % GO TRHOUGH DATA TO SEE EXCEPTIONS (use Data Cursor)
        figure; 
        c = 1;
        for cid = 1:3
            for speedid = 1:3
                subplot(3,3,c);
                imagesc(S(cid,speedid).dx_lock_int',[-10 50]); colorbar;
                title(sprintf('COND %d, SPEED %d\n', cid, speedid));
                c = c+1;
            end
        end
        
        save(datafile,'S'); % too much data, compression
        fprintf('go to next sub, save %s\n', datafile);
        
    end
    
end 
                    
%% if want to rerun the analysis
StartSample = 1;
figure(1);
    for subcnt = 1:length(sub_list)

%keyboard
        datafile = [int2str( sub_list(subcnt) ) '.mat'];
        load(datafile,'S'); % too much data, compression

        % Check if field exists, otherwise load data
        if ~isfield(S, 'PHASEx')

            for cid = 1:3
                for speedid = 1:3            

                    % not on the locked traces but the full trace 
                    [ntrl,~] = size(S(cid,speedid).ttdx);            

                    nF = 100; % n frequencies stored

                    % target: for every different "cycle_ms"
                    EndSample = length(S(cid,speedid).ttdx(1,:))-100;

                    lengthz = length(S(cid,speedid).ttdx(1,StartSample:EndSample));

                    % ONLY STORE THE TARGET FREQUENCY IN S
                    % USE TMP struct for the grand average
                    S(cid,speedid).PHASEx = NaN.*ones(1,lengthz, ntrl);
                    S(cid,speedid).COHx = NaN.*ones(1,lengthz, ntrl);

                    S(cid,speedid).PHASEy = NaN.*ones(1,lengthz, ntrl);
                    S(cid,speedid).COHy = NaN.*ones(1,lengthz, ntrl);

                    % FULL MATRIX
                    TMP(cid,speedid).PHASEx = NaN.*ones(nF,lengthz, ntrl);
                    TMP(cid,speedid).COHx = NaN.*ones(nF,lengthz, ntrl);

                    TMP(cid,speedid).PHASEy = NaN.*ones(nF,lengthz, ntrl);
                    TMP(cid,speedid).COHy = NaN.*ones(nF,lengthz, ntrl);

                    TMP(cid,speedid).dx = NaN.*ones(1,lengthz, ntrl);
                    TMP(cid,speedid).dy = NaN.*ones(1,lengthz, ntrl);

                    % target
                    TMP(cid,speedid).tx = NaN.*ones(1,lengthz, ntrl);
                    TMP(cid,speedid).ty = NaN.*ones(1,lengthz, ntrl);

                    for n = 1:ntrl
                        
                        tx = S(cid,speedid).ttdx(n,StartSample:EndSample)';
                        ty = S(cid,speedid).ttdy(n,StartSample:EndSample)';

                        x = S(cid,speedid).ointdx(n,StartSample:EndSample); %StartSample:(StartSample+lengthz-1));
                        y = S(cid,speedid).ointdy(n,StartSample:EndSample); %StartSample:(StartSample+lengthz-1));
                        t = StartSample:(lengthz+StartSample-1);
                        
                        % find the first NaN and cut there ...
                        % pad for saving
                        x(1) = 0;
                        y(1) = 0;

                        % FIND THE LAST Non NaN
                        l = find(~isnan(x)); 
                        ln = l(end);
                        
                        % First non-NaN                       
                        ls = l(1);
                        
                        % INTERPOLATION NECESSARY TO BRIDGE THE GAP BETWEEN
                        % BLINKS
                        if sum(isnan(x(ls:ln))) || sum(isnan(y(ls:ln)))
                              % INTERPOLATION OF NaNs
                              xnonan =inpaint_nans(x(ls:ln),5); % method 5: average neighbour
                              ynonan =inpaint_nans(y(ls:ln),5);
                              %clf; 
                              %hold on; plot(xnonan,'b'); hold on; plot(x,'r'); 
%                              hold on; plot(ynonan,'b'); hold on; plot(y,'r'); 

                        else
                            xnonan = x(ls:ln);
                            ynonan = y(ls:ln);
                        end
                        try

                            [WCOHx,WCSx,F,COI] = fCoherenceAnalysis(xnonan, tx(ls:ln), 9, 32);
                            [WCOHy,WCSy,F,COI] = fCoherenceAnalysis(ynonan, ty(ls:ln), 9, 32);
                            
                            % save ONLY THE FREQUENCY OF THE TARGET:
                            d = abs(F-TRAJ(cid,speedid).PK_X);
                            FtoSave = find(d==min(d));

                            % save signal for checking and plots: x y
                            TMP(cid,speedid).dx(:,ls:ln,n) = x(1:ln); 
                            TMP(cid,speedid).dy(:,ls:ln,n) = y(1:ln); 

                            % target
                            TMP(cid,speedid).tx(:,ls:ln,n) = tx(1:ln);
                            TMP(cid,speedid).ty(:,ls:ln,n) = ty(1:ln);

                            TMP(cid,speedid).PHASEx(:,ls:ln,n) = -angle(WCSx(end-nF+1:end,:));
                            TMP(cid,speedid).COHx(:,ls:ln,n) = WCOHx(end-nF+1:end,:);

                            TMP(cid,speedid).PHASEy(:,ls:ln,n) = -angle(WCSy(end-nF+1:end,:));
                            TMP(cid,speedid).COHy(:,ls:ln,n) = WCOHy(end-nF+1:end,:);

                            S(cid,speedid).PHASEx(1,ls:ln,n) = -angle(WCSx(FtoSave,:))';
                            S(cid,speedid).COHx(1,ls:ln,n) = WCOHx(FtoSave,:);

                            S(cid,speedid).PHASEy(1,ls:ln,n) = -angle(WCSy(FtoSave,:));
                            S(cid,speedid).COHy(1,ls:ln,n) = WCOHy(FtoSave,:);

                            TMP(cid,speedid).t = t;
                            TMP(cid,speedid).F = F(end-nF+1:end);

                            TMP(cid,speedid).COI = COI;

                            S(cid,speedid).t = t;
                            S(cid,speedid).F = F(FtoSave);

                            S(cid,speedid).COI = COI;

                            %keyboard
                            fprintf('trial %d cid %d speedid %d subcnt %d\n', n, cid, speedid, subcnt);   

                        catch me

                            me
                            %keyboard

                        end

                    end
                end
            end

            % AVERAGE MAP
            t = 1;
            for c = 1:3
                for s = 1:3

                    % CAN ALSO GIVE CONFIDENCE LIMITS
                    GPHASEx = circ_mean(TMP(c,s).PHASEx, [], 3); % CIRCULAR DATA! CIRC_MEAN!
                    GPHASEy = circ_mean(TMP(c,s).PHASEy, [], 3); % CIRCULAR DATA! CIRC_MEAN!
                    GCOHx = nanmean(TMP(c,s).COHx(:,:,:), 3);
                    GCOHy = nanmean(TMP(c,s).COHy(:,:,:), 3);            

                    S(c,s).GPHASEx = GPHASEx; %circ_mean(S(c,s).PHASEx, [], 3); % CIRCULAR DATA! CIRC_MEAN!
                    S(c,s).GPHASEy = GPHASEy; %circ_mean(S(c,s).PHASEy, [], 3); % CIRCULAR DATA! CIRC_MEAN!
                    S(c,s).GCOHx = GCOHx; %nanmean(S(c,s).COHx(:,:,:), 3);
                    S(c,s).GCOHy =  GCOHy; %nanmean(S(c,s).COHx(:,:,:), 3);            

                    t = t+1;

                end
            end    
            save(datafile,'S');
        end
end

%% GROUP ACCROSS SUBJECTS
G=struct([]);
subn = 1;

for sub = sub_list %subn = 1;
    
    load([int2str(sub) '.mat']);

    for c = 1:3
        for s = 1:3

            fprintf('subject %d, type %d, speed %d\n', subn, c, s);

            % TBD: field names
            G(c,s).GPHASEx(:,:,subn)= S(c,s).GPHASEx;
            G(c,s).GPHASEy(:,:,subn)= S(c,s).GPHASEy;
            G(c,s).GCOHx(:,:,subn)= S(c,s).GCOHx;        
            G(c,s).GCOHy(:,:,subn)= S(c,s).GCOHy;        

            % accross trials (to GENERIC FCT)
            % Get size of struct
            
            % Init max size
            maxt = 4000;
            maxtrl = 50;
            G(c,s).PHASEx(:,:,subn)= NaN.*ones(maxtrl,maxt);
            G(c,s).PHASEy(:,:,subn)= NaN.*ones(maxtrl,maxt);
            G(c,s).COHx(:,:,subn)= NaN.*ones(maxtrl,maxt);
            G(c,s).COHy(:,:,subn)=  NaN.*ones(maxtrl,maxt);       

            [v b m ] = size(S(c,s).PHASEx);
            % recode to have trials x samples
            % shiftdim: shift the leading singleton dimension, or squeeze to the same effect...
            G(c,s).PHASEx(1:m,1:b,subn)= squeeze(S(c,s).PHASEx)';             
            G(c,s).PHASEy(1:m,1:b,subn)= squeeze(S(c,s).PHASEy)';        
            G(c,s).COHx(1:m,1:b,subn)= squeeze(S(c,s).COHx)';        
            G(c,s).COHy(1:m,1:b,subn)= squeeze(S(c,s).COHy)';        

        end
    end
    subn = subn + 1;
end

%% PLOT THE AVERAGE COHERENCE AND PHASE PLOT
% Trial by trial effects ...

[WCOHx,WCSx,F,COI] = fCoherenceAnalysis(TRAJ(1,1).dx, TRAJ(1,1).dx, 9, 32);

figure;
%DSFigDef(2,[1. 1. 1],1);

t = 1;
range = 100;
Freq = F(end-range:end);

for c = 1:3
    for s = 1:3
        subplot(3,3,t);
        imagesc(circ_mean(G(c,s).GPHASEx,[],3),[-pi pi]); colorbar;
        title(sprintf('Phase X: c %d s %d, F %2.2f\n', c,s, S(c,s).F));
        hold on; plot(range-(S(c,s).COI.*(range/Freq(1))),'--w');
        
        % target principal freq:         
        d = abs(Freq-S(c,s).F);
        FtoSave = find(d==min(d));

        set(gca, 'YTick',1:10:range,'YTickLabel',fix(Freq(1:10:range)*10)/10);
        hold on; plot([0 4000],[FtoSave FtoSave],'--r','LineWidth',1);

        t = t+1;
    end
end

figure; 
%DSFigDef(2,[1. 1. 1],1);
t = 1;
for c = 1:3
    for s = 1:3
        subplot(3,3,t);
        imagesc(circ_mean(G(c,s).GPHASEy,[],3),[-pi pi]); colorbar;
        title(sprintf('Phase Y: c %d s %d, F %2.2f\n', c,s, S(c,s).F));
        hold on; plot(range-(S(c,s).COI.*(range/Freq(1))),'--w');
        
        % target principal freq:         
        d = abs(Freq-S(c,s).F);
        FtoSave = find(d==min(d));

        set(gca, 'YTick',1:10:range,'YTickLabel',fix(Freq(1:10:range)*10)/10);
        hold on; plot([0 4000],[FtoSave FtoSave],'--r','LineWidth',1);

        t = t+1;
    end
end

figure;
%DSFigDef(2,[1. 1. 1],1);
t = 1;
for c = 1:3
    for s = 1:3
        subplot(3,3,t);
        imagesc(nanmean(G(c,s).GCOHx,3),[0 1]); colorbar;

        title(sprintf('COH X: c %d s %d, F %2.2f\n', c,s, S(c,s).F));
        hold on; plot(range-(S(c,s).COI.*(range/Freq(1))),'--w');
        
        % target principal freq:         
        d = abs(Freq-S(c,s).F);
        FtoSave = find(d==min(d));

        set(gca, 'YTick',1:10:range,'YTickLabel',fix(Freq(1:10:range)*10)/10);
        hold on; plot([0 4000],[FtoSave FtoSave],'--r','LineWidth',1);

        t = t+1;
    end
end
figname = 'Fig3a';
print(gcf, '-depsc2', [pwd '\..\..\figs\' figname '.eps']);
print(gcf, '-dpdf',  [pwd '\..\..\figs\' figname '.pdf']);     
%plot2svg([pwd '\..\..\figs\' figname '.svg']);

figure;
%DSFigDef(2,[1. 1. 1],1);
t = 1;
for c = 1:3
    for s = 1:3
        subplot(3,3,t);
        imagesc(nanmean(G(c,s).GCOHy,3),[0 1]); colorbar;
               
        title(sprintf('COH Y: c %d s %d, F %2.2f\n', c,s, S(c,s).F));
        hold on; plot(range-(S(c,s).COI.*(range/Freq(1))),'--w');
        
        % target principal freq:         
        d = abs(Freq-S(c,s).F);
        FtoSave = find(d==min(d));

        set(gca, 'YTick',1:10:range+1,'YTickLabel',fix(Freq(1:10:range+1)*10)/10);
        hold on; plot([0 4000],[FtoSave FtoSave],'--r','LineWidth',1);
        
        t = t+1;
    end
end

%wcoherence(nanmean(G(c,s).GCOHy,3),1:3000,F,coi,'Seconds','Hz');

figname = 'Fig3b';
print(gcf, '-depsc2', [pwd '\..\..\figs\' figname '.eps']);
print(gcf, '-dpdf',  [pwd '\..\..\figs\' figname '.pdf']);     
%plot2svg([pwd '\..\..\figs\' figname '.svg']);

%% PLOT: (new fct, go through dims to make subplots ...)
% Trial by trial effects ...
figure;
t = 1;
for c = 1:3
    for s = 1:3
        subplot(3,3,t);
        imagesc(circ_mean(G(c,s).PHASEx,[],3),[-pi pi]); colorbar;
        title(sprintf('Phase X: c %d s %d\n', c,s));

        [szx ~] = size(TRAJ(c,s).dx);
        
        xlim([0 szx]);
        
        t = t+1;
    end
end

figure;
t = 1;
for c = 1:3
    for s = 1:3
        subplot(3,3,t);
        imagesc(circ_mean(G(c,s).PHASEy,[],3),[-pi pi]); colorbar;
        title(sprintf('Phase Y: c %d s %d\n', c,s));
        [szx ~] = size(TRAJ(c,s).dx);
        
        xlim([0 szx]);

        t = t+1;
    end
end

%Freq=  F((end-nF):end);
     
figure;
t = 1;
for c = 1:3
    for s = 1:3
        subplot(3,3,t);
        imagesc(nanmean(G(c,s).COHx,3),[0 1]); colorbar;
  
        title(sprintf('COH X: c: %d, s: %d\n', c,s));
        %axis([0 3000 0.2 .6]);
        [szx ~] = size(TRAJ(c,s).dx);
        
        xlim([0 szx]);
      
        t = t+1;
    end
end

figure;
t = 1;
for c = 1:3
    for s = 1:3
        subplot(3,3,t);
        imagesc(nanmean(G(c,s).COHy,3),[0 1]); colorbar;
        
        title(sprintf('COH Y: c: %d, s: %d\n', c,s));
        %%axis([0 3000 0.2 .6]);
        [szx ~] = size(TRAJ(c,s).dx);
        
        xlim([0 szx]);
        
        t = t+1;
    end
end

%% FOR X and Y show average coherence accross trials ...
figure;
p = panel();
p.pack(2, 3);
p.de.margin = 5;
p.fontsize = 10;
p.fontname = 'Arial';

% set margins
%p(1,1).marginbottom = 12;
%p(2).marginleft = 20;

% and some properties

cols = 'rgb';
%t = 1;
for s = 1:3
    for c = 1:3
        %subplot(2,3,t);
        p(1,s).select();
        p(1,s).margin = [0 0 5 0];
        if s == 1
            p(1,s).ylabel('Horizontal Coherence');
            set(gca, 'xtick', []);

        else
            set(gca, 'xtick', [], 'ytick', []);
        end
        
        % average trials
        x = nanmean(G(c,s).COHx(:,:,:), 2);
        % then average across subs

        hold on; errorbar(nanmean(x,3), stderr(x,3),'-','color', COLMAP(c,:));
        
        % smoothed version: 
        y = nanmean(x,3);
        tx = (1:length(y))';
        id = ~isnan(y);
        yy2 = smooth(tx(id),y(id),.4,'lowess');
        hold on; plot(yy2,'color',COLMAP(c,:),'LineWidth',2);
        box;
        ylim([.2 .9]);
    end
    text(5,.85, sprintf('%2.2f deg/sec %d\n', mean(TRAJ(1,s).dx)));

end

t = 1;
for s = 1:3
    for c = 1:3
        %subplot(2,3,t+3);
        p(2,s).select();
        p(2,s).margin = [0 0 5 0];
        p(2,s).xlabel('Trial number');
        if s == 1
            p(2,s).ylabel('Vertical Coherence');
            %set(gca, 'xtick', [], 'ytick', []);
        else
            set(gca, 'ytick', []);
        end
        
        x = nanmean(G(c,s).COHy(:,:,:),2);
        
        % then average across subs
        hold on; errorbar(nanmean(x,3),stderr(x,3),'-','color', COLMAP(c,:));
        box;
        y = nanmean(x,3);
        tx = (1:length(y))';
        id = ~isnan(y);
        yy2 = smooth(tx(id),y(id),.4,'lowess');
        hold on; plot(yy2,'color', COLMAP(c,:),'LineWidth',2);
        ylim([.2 .9]);
        
    end
    if t==1
       % legend('SouthWest', {'upright','inverted','scrambled'});
    end
    text(5,.85, sprintf('%2.2f deg/sec %d\n', mean(TRAJ(1,s).dx)));
    t=t+1;
end

figname = 'Fig4';
print(gcf, '-depsc2', [pwd '\..\..\figs\' figname '.eps']);
print(gcf, '-dpdf',  [pwd '\..\..\figs\' figname '.pdf']);     
%plot2svg([pwd '\..\..\figs\' figname '.svg']);

%% PHASE
figure; 
%DSFigDef(2,[1. 1. 1],1);

% TO CHECK IT IS ACCURATE: PROBLEM, HAVE ZERO'S AT TIMES INSTEAD OF NaN
%INTERVAL = 1000:2000;
t = 1;
cols = 'rgb';
figure;
for s = 1:3
    for c = 1:3
        subplot(2,3,t);
        % average trials
        %circ_std(x,[],[],3)./sqrt(21)

        % average along time intervals, then subs
        x = circ_mean(G(c,s).PHASEx(1:40,:,:), [], 2);
        
        % assign NaN if no data within the trial
%        x(x==0) = NaN;
        
        % then average across subs
        subx= circ_mean(x,[],3); % confidence interval
        hold on; errorbar(subx, circ_std(x,[],[],3)./sqrt(21),'-', 'color', COLMAP(c,:));
        set(gca,'YTick',[-pi -pi/2 0 pi/2 pi],'YTickLabel',{'-pi' '-pi/2' '0' 'pi/2' 'pi'});
        title(sprintf('PHASE X: c: %d, s: %d\n', c,s));
        
        hold on; plot([0 40],[0 0],'--k','LineWidth',2);

%         y = circ_mean(x,[],3);
%         tx = (1:length(y))';
%         id = ~isnan(y);
%         yy2 = smooth(tx(id),y(id),.4,'lowess');
%         hold on; plot(yy2,cols(c),'LineWidth',2);
        axis([0 40 -pi pi]);
    end
    t=t+1;
end

figure

t = 1;
for s = 1:3
    for c = 1:3
        subplot(2,3,t+3);
        x = circ_mean(G(c,s).PHASEy(1:40,:,:), [], 2);
        % then average across subs
        % assign NaN if no data within the trial
        x(x==0) = NaN;

        % recode: everything that is positive, recode as x = -x - pi.
%        posid = x > 0;
%        x(posid)= (-pi+x(posid)) -pi;
        hold on; errorbar(circ_mean(x,[],3), circ_std(x,[], [], 3)./sqrt(21),'-', 'color', COLMAP(c,:));
        title(sprintf('PHASE y: c: %d, s: %d\n', c,s));
        
%         y = circ_mean(x,[],3);
%         tx = (1:length(y))';
%         id = ~isnan(y);
%         yy2 = smooth(tx(id),y(id),.4,'lowess');
%         hold on; plot(yy2,cols(c),'LineWidth',2);
%         hold on; plot([0 40],[0 0],'--k','LineWidth',2);

% relabel
%[-pi -pi/2 0 pi/2 pi]+pi
        set(gca,'YTick',[-pi -pi/2 0 pi/2 pi],'YTickLabel',{'0' '-pi/2' '-pi' 'pi' 'pi/2'});
        axis([0 40 -pi pi]);
    end
    t = t+1;    
end

%% Plot average angle per trial, quiver(x,y,u,v)
% FIGURE 5

figure;
p = panel();
p.pack(2, 3); 
p.de.margin = 1;
p.fontsize = 10;
p.fontname = 'Arial';

% Show different lags accross time for X and Y
t = 1;
for s = 1:3
    p(1,s).pack(3,1);
    for c = 1:3
        %subplot(2,3,t);
        % average trials
        %circ_std(x,[],[],3)./sqrt(21)
        p(1,s,c,1).select();
        
        % average along time intervals, then subs
        x = circ_mean(G(c,s).PHASEx(1:40,:,:), [], 2);
        
        x(x==0) = NaN;

        % then average across subs
        subx= circ_mean(x,[],3); % confidence interval
        angle = subx;
        % determine the points to place the arrows, arbitrary here trial numbers
        [x,y] = meshgrid(1:40,0);
        norm = 2;
        u = cos(angle).*norm;
        v = sin(angle).*norm;

        quiver(x,y,u',v', 'color', COLMAP(c,:));
        axis([0 40 -10 10]);
        set(gca,'ytick',[],'xtick',[]);
        
%        axis([0 40 -pi pi]);
    end
    t=t+1;
end

t = 1;
for s = 1:3
    p(2,s).pack(3,1);
    for c = 1:3
        p(2,s,c,1).select();
        
        %subplot(2,3,t+3);
        x = circ_mean(G(c,s).PHASEy(1:40,:,:), [], 2);
        % then average across subs
        % assign NaN if no data within the trial
        x(x==0) = NaN;
        subx= circ_mean(x,[],3); % confidence interval

        angle = subx;
        % determine the points to place the arrows, arbitrary here trial numbers
        [x,y] = meshgrid(1:40,0);
        norm = 2;
        u = cos(angle).*norm;
        v = sin(angle).*norm;

        quiver(x,y,u',v','color', COLMAP(c,:));
        axis([0 40 -10 10]);

        set(gca,'ytick',[],'xtick',[]);
        
    end
    t = t+1;    
end
figname = 'Fig5';
print(gcf, '-depsc2', [pwd '\..\..\figs\' figname '.eps']);
print(gcf, '-dpdf',  [pwd '\..\..\figs\' figname '.pdf']);     

%% Phase shift accross the trial (quiver plot)
% FIGURE 6
figure;
p = panel();
p.pack(2, 3); 
p.de.margin = 1;
p.fontsize = 10;
p.fontname = 'Arial';

% Show different lags accross time for X and Y
t = 1;
for s = 1:3
    p(1,s).pack(3,1);
    for c = 1:3
        %subplot(2,3,t);
        % average trials
        %circ_std(x,[],[],3)./sqrt(21)
        p(1,s,c,1).select();
        
        % average along time intervals, then subs
        x = circ_mean(G(c,s).PHASEx(:,:,:), [], 1);
        
        x(x==0) = NaN;

        % then average across subs
        subx= circ_mean(x,[],3); % confidence interval
        angle = subx;
        
        % 20 bins
        bin_angle = average_time_bins(angle(1:length(TRAJ(1,s).dx)),20);
        
        % determine the points to place the arrows, arbitrary here trial numbers
        [x,y] = meshgrid(1:20,0);
        norm = 2;
        u = cos(bin_angle).*norm;
        v = sin(bin_angle).*norm;

        quiver(x,y,u,v, 'color', COLMAP(c,:));
        ylim([-5 5]);
%        keyboard
        set(gca,'ytick',[],'xtick',(1:4:20),'xtickLabel',[]);
    end
    t=t+1;
end

t = 1;
for s = 1:3
    p(2,s).pack(3,1);
    for c = 1:3
        %subplot(2,3,t);
        % average trials
        %circ_std(x,[],[],3)./sqrt(21)
        p(2,s,c,1).select();
        
        % average along time intervals, then subs
        x = circ_mean(G(c,s).PHASEy(:,:,:), [], 1);
        
        x(x==0) = NaN;

        % then average across subs
        subx= circ_mean(x,[],3); % confidence interval
        angle = subx;
        
        % 20 bins
        bin_angle = average_time_bins(angle(1:length(TRAJ(1,s).dx)),20);
        
        % determine the points to place the arrows, arbitrary here trial numbers
        [x,y] = meshgrid(1:20,0);
        norm = 2;
        u = cos(bin_angle).*norm;
        v = sin(bin_angle).*norm;

        quiver(x,y,u,v, 'color', COLMAP(c,:));
        ylim([-5 5]);
%        keyboard
        set(gca,'ytick',[],'xtick',(1:4:20),'xtickLabel',(1:4:20)*4000/20);
    end
    t = t+1;    
end

figname = 'Fig6';
print(gcf, '-depsc2', [pwd '\..\..\figs\' figname '.eps']);
print(gcf, '-dpdf',  [pwd '\..\..\figs\' figname '.pdf']);     

%% POLAR PLOT
% Figure 7
figure;
p = panel();
p.pack(2, 3); 
p.de.margin = 1;
p.fontsize = 10;
p.fontname = 'Arial';
%title(sprintf('Horizontal eye movement, %s deg/sec target',Ctrl.vels(s)));

t = 1;
for s = 1:3
    
   % subplot(2,3,t);
%    polarhistogram([],30,'FaceColor',cols(c));
   % x = circ_mean(G(1,s).PHASEy(1:40,INTERVAL,:), [], 2);

    %circ_plot(circ_mean(x,[],3),'pretty','ro',true);
    p(1,s).select();
    p(1,s).margin = [10 10 10 10]; % right, xx, left,  
    for c = 1:3
        x = circ_mean(G(c,s).PHASEx(1:40,:,:), [], 2);
        % then average across subs
        % assign NaN if no data within the trial

        hold on; circ_plot(circ_mean(x,[],3),'pretty','ro',true,'color',COLMAP(c,:)),
        set(gca, 'Visible','off')
    end
    title(sprintf('Horizontal EM, %2.2f deg/sec %d\n', mean(TRAJ(1,s).dx)));

    t = t+1;    
end

t = 1;
for s = 1:3
    
    %subplot(2,3,t+3);
    p(2,s).select();
    p(2,s).margin = [10 10 10 10]; % right, xx, left,      
    for c = 1:3
        x = circ_mean(G(c,s).PHASEy(1:40,:,:), [], 2);

        hold on; circ_plot(circ_mean(x,[],3),'pretty','ro',true,'color',COLMAP(c,:)),
        set(gca, 'Visible','off')
    
    end
    title(sprintf('Vertical EM, %2.2f deg/sec %d\n', mean(TRAJ(1,s).dx)));

    t = t+1;    
end

figname = 'Fig7';
print(gcf, '-depsc2', [pwd '\..\..\figs\' figname '.eps']);
print(gcf, '-dpdf',  [pwd '\..\..\figs\' figname '.pdf']);     

%% NEED TO SHOW IN THOSE PLOTS THE POINTS FROM THE CONE OF CONFIDENCE(?)
t = 1;
cols = 'rgb';
figure;
for s = 1:3
    for c = 1:3
        subplot(2,3,t);
        % average trials
        x = nanmean(G(c,s).COHx, 1);
        % then average across subs
        hold on; errorbar(nanmean(x, 3), stderr(x,3),['-' cols(c)]);
        title(sprintf('COH X: c: %d, s: %d\n', c,s));
        axis([500 2000 0 1]);
    end
    t=t+1;
end
t = 1;
for s = 1:3
    for c = 1:3
        subplot(2,3,t+3);
        x = nanmean(G(c,s).COHy, 1);
        
        % then average across subs
        hold on; errorbar(mean(x,3),stderr(x,3),['-' cols(c)]);
        title(sprintf('COH Y: c: %d, s: %d\n', c,s));
        axis([500 2000 0 1]);
    end
    t = t+1;    
end
%%
t = 1;
cols = 'rgb';
figure;
for s = 1:3
    for c = 1:3
        subplot(2,3,t);
        % average trials
        x = circ_mean(G(c,s).PHASEx, [], 1);
        x(x==0) = NaN;

        hold on; errorbar(circ_mean(x,[], 3), circ_std(x,[],[],3)./sqrt(21),['-' cols(c)]);
        title(sprintf('Phase X\n', c,s));
        axis([0 3000 -pi pi]);
    end
    t = t + 1;
end
t = 1;
for s = 1:3
    for c = 1:3
        subplot(2,3,t+3);
        x = circ_mean(G(c,s).PHASEy, [], 1);
        x(x==0) = NaN;

        % then average across subs
        hold on; errorbar(circ_mean(x,[],3),circ_std(x,[],[],3)./sqrt(21),['-' cols(c)]);
        
        % Express all as positive values:
        title(sprintf('Phase Y %d %d\n', c,s));
        axis([0 3000 -pi pi]);
        
    end
    t = t+1;    
end

%% FIGURE 8: Perception
% 
t = 1;
for n = sub_list

    Trl.sub = n;
    Man.log = sprintf('%se%ds%d/e%ds%d.log', Man.datpath, Trl.exp, Trl.sub, Trl.exp, Trl.sub);
    DATA = textread(Man.log);
    MDATA = DATA(DATA(:,Col.train)==0 & DATA(:,Col.task)==0,:); 
    sPercMatches;
    
    datafile = [int2str(n) '.mat'];
    load(datafile,'S'); % too much data, compression
    fprintf('%s\n', datafile);
    % for anova and lineplot
    for c=1:3
        for s=1:3
            
            M(t,:) = [abs(MATCH(c,s).m) c s n];
            GAIN_INT=500:length(TRAJ(c,s).dx);
            GG(t,:) = [nanmean(nanmean(S(c,s).ointdx(:,GAIN_INT),2))./mean(TRAJ(c,s).dx) c s n];
%            G(t,:) = [nanmean(S(c,s).gainx) c s n];
            t = t + 1;
        end
    end    
end

%% FIGURE 8, matches and FIGURE 9 gain
cc.match=1;
cc.type=2;
cc.speed=3;
cc.sub=4;
title('Perceptual matches');
id = ~isnan(M(:,cc.match));
barplot(M(id,cc.type), M(id,cc.match), 'split', M(id,cc.speed),'style_rainbow','leg',{'low','med','fast'});%,'linecolor', COLMAP);

xlabel('Condition');
ylim([0.5 1.1]);

figname = 'Fig8';
print(gcf, '-depsc2', [pwd '\..\..\figs\' figname '.eps']);
print(gcf, '-dpdf',  [pwd '\..\..\figs\' figname '.pdf']);     

% Save for Kyle to do an ANOVA, follow-up t-tests. Write-up.
dlmwrite('anova_matches.txt',M(id,:),'delimiter',';');
% TBD: Correlate perception and coherence estimates ...

% FIGURE 9: Gain
id = ~isnan(GG(:,1));
barplot(GG(id,cc.type), GG(id,1), 'split', M(id,cc.speed),'style_rainbow','leg',{'low','med','fast'});%,'linecolor', COLMAP);

ylim([.9 1.0]);

figname = 'Fig9';
print(gcf, '-depsc2', [pwd '\..\..\figs\' figname '.eps']);
print(gcf, '-dpdf',  [pwd '\..\..\figs\' figname '.pdf']);     

% Save for Kyle to do an ANOVA, follow-up t-tests. Write-up.
dlmwrite('gain.txt',M(id,:),'delimiter',';');

%% FIGURE 1: show perceptual task possible matches ...
figure;

% Figure 1a: plot dots at different moments in time
title('bounciness manipulation');
for n = 0.25:0.25:1.5

    Trl.match = n;
    Trl.cond = 1;
    Trl.speed = Ctrl.speed(speed);
    Trl.cycle0 = 0; % zero shift! ; %cycle0 = Trl.cycle0; %DATA(Trl.trl,Col.cycle0);

    Trl.dir = 0;

    % assign traj structure (Trl.cond)
    traj = fGetWalkerTraj(Trl,Col,Ctrl);  

    hold on; plot((1:length(traj.dy))./1000,traj.dy,'color',1-[.5 .5 .5].*n,'LineWidth',3);
end

legend({'0.25','0.5','0.75','1','1.25','1.5'}); %:0.25:1.5);
hold on;
%Trl.match = 1;

%traj = fGetWalkerTraj(Trl,Col,Ctrl);
%hold on; plot((1:length(traj.dy))./1000,traj.dy,'r','LineWidth',2);
xlabel('time [s]');
ylabel('Vertical velocity [deg/sec]');
axis([0 2 -4 4]);

figname = 'Fig2b';
print(gcf, '-depsc2', [pwd '\..\..\figs\' figname '.eps']);
print(gcf, '-dpdf',  [pwd '\..\..\figs\' figname '.pdf']);     

% Figure 1a: show the dots traj accross the screen... X, Y
% plot(X,Y); plot(hipx, hipy);

% TBD: NEED TO BE ABLE TO LOOK AT THE GAIT CYCLE IN RELATION TO THE HIP
% VELOCITY...
%%
traj = fGetWalkerTraj(Trl,Col,Ctrl);
figure; 
hip_dot = 13;

figure


for t = 1:20:1000

    % connect dots...
    
%    hold on; plot(traj.allx(t,1), traj.ally(t,1),'-ok','MarkerSize',2);
    hold on; plot(traj.allx(t,[1:2 9:10]), traj.ally(t,[1:2 9:10]),'-ok','MarkerFaceColor','k','MarkerSize',2);
    hold on; plot(traj.allx(t,3:5), traj.ally(t,3:5),'-ok','MarkerFaceColor','k','MarkerSize',2);
    hold on; plot(traj.allx(t,6:8), traj.ally(t,6:8),'-ok','MarkerFaceColor','k','MarkerSize',2);
    hold on; plot(traj.allx(t,13:15), traj.ally(t,13:15),'-og','MarkerFaceColor','g','MarkerSize',2);
    hold on; plot(traj.allx(t,[13 11 12]), traj.ally(t,[13 11 12]),'-ok','MarkerFaceColor','k','MarkerSize',2);
    hold on; plot(traj.allx(t,13), traj.ally(t,13),'-or','MarkerFaceColor','r','MarkerSize',2);
    
end
hold on; plot(traj.allx(1000:end,hip_dot), traj.ally(1000:end,hip_dot),'-r','LineWidth',2);
xlabel('Horizontal Position [deg]');
ylabel('Vertical Position [deg]');

pbaspect([1.5 1 1])

figname = 'Fig1a';
print(gcf, '-depsc2', [pwd '\..\..\figs\' figname '.eps']);
print(gcf, '-dpdf',  [pwd '\..\..\figs\' figname '.pdf']);     

%% Figure 1d-e: single trace examples

% SHOW THE PHASE INFO ...

clear angle;
datafile = [int2str( sub_list(2) ) '.mat'];
load(datafile,'S'); % too much data, compression

cid = 1;
speedid =1;
n = 6;
StartSample= 1;

[ntrl,~] = size(S(cid,speedid).ttdx);            
nF = 100; % n frequencies stored
EndSample = length(S(cid,speedid).ttdx(1,:))-100;
lengthz = length(S(cid,speedid).ttdx(1,StartSample:EndSample));

S(cid,speedid).PHASEx = NaN.*ones(1,lengthz, ntrl);
S(cid,speedid).COHx = NaN.*ones(1,lengthz, ntrl);
S(cid,speedid).PHASEy = NaN.*ones(1,lengthz, ntrl);
S(cid,speedid).COHy = NaN.*ones(1,lengthz, ntrl);

% FULL MATRIX
TMP(cid,speedid).PHASEx = NaN.*ones(nF,lengthz, ntrl);
TMP(cid,speedid).COHx = NaN.*ones(nF,lengthz, ntrl);

TMP(cid,speedid).PHASEy = NaN.*ones(nF,lengthz, ntrl);
TMP(cid,speedid).COHy = NaN.*ones(nF,lengthz, ntrl);

TMP(cid,speedid).dx = NaN.*ones(1,lengthz, ntrl);
TMP(cid,speedid).dy = NaN.*ones(1,lengthz, ntrl);

% target
TMP(cid,speedid).tx = NaN.*ones(1,lengthz, ntrl);
TMP(cid,speedid).ty = NaN.*ones(1,lengthz, ntrl);

tx = S(cid,speedid).ttdx(n,StartSample:EndSample)';
ty = S(cid,speedid).ttdy(n,StartSample:EndSample)';

x = S(cid,speedid).ointdx(n,StartSample:EndSample); %StartSample:(StartSample+lengthz-1));
y = S(cid,speedid).ointdy(n,StartSample:EndSample); %StartSample:(StartSample+lengthz-1));
t = StartSample:(lengthz+StartSample-1);
                        
% find the first NaN and cut there ...
% pad for saving
x(1) = 0;
y(1) = 0;

% FIND THE LAST Non NaN
l = find(~isnan(x)); 
ln = l(end);
                        
% First non-NaN                       
ls = l(1);
                        
% INTERPOLATION NECESSARY TO BRIDGE THE GAP BETWEEN
% BLINKS
if sum(isnan(x(ls:ln))) || sum(isnan(y(ls:ln)))
      % INTERPOLATION OF NaNs
      xnonan =inpaint_nans(x(ls:ln),5); % method 5: average neighbour
      ynonan =inpaint_nans(y(ls:ln),5);
      %clf; 
      %hold on; plot(xnonan,'b'); hold on; plot(x,'r'); 
%                              hold on; plot(ynonan,'b'); hold on; plot(y,'r'); 

else
    xnonan = x(ls:ln);
    ynonan = y(ls:ln);
end

[WCOHx,WCSx,F,COI] = fCoherenceAnalysis(xnonan, tx(ls:ln), 9, 32);
[WCOHy,WCSy,F,COI] = fCoherenceAnalysis(ynonan, ty(ls:ln), 9, 32);

% save ONLY THE FREQUENCY OF THE TARGET:
d = abs(F-TRAJ(cid,speedid).PK_X);
FtoSave = find(d==min(d));

% save signal for checking and plots: x y
TMP(cid,speedid).dx(:,ls:ln,n) = x(1:ln); 
TMP(cid,speedid).dy(:,ls:ln,n) = y(1:ln); 

% target
TMP(cid,speedid).tx(:,ls:ln,n) = tx(1:ln);
TMP(cid,speedid).ty(:,ls:ln,n) = ty(1:ln);

TMP(cid,speedid).PHASEx(:,ls:ln,n) = -angle(WCSx(end-nF+1:end,:));
TMP(cid,speedid).COHx(:,ls:ln,n) = WCOHx(end-nF+1:end,:);

TMP(cid,speedid).PHASEy(:,ls:ln,n) = -angle(WCSy(end-nF+1:end,:));
TMP(cid,speedid).COHy(:,ls:ln,n) = WCOHy(end-nF+1:end,:);

S(cid,speedid).PHASEx(1,ls:ln,n) = -angle(WCSx(FtoSave,:))';
S(cid,speedid).COHx(1,ls:ln,n) = WCOHx(FtoSave,:);

S(cid,speedid).PHASEy(1,ls:ln,n) = -angle(WCSy(FtoSave,:));
S(cid,speedid).COHy(1,ls:ln,n) = WCOHy(FtoSave,:);

TMP(cid,speedid).t = t;
TMP(cid,speedid).F = F(end-nF+1:end);

TMP(cid,speedid).COI = COI;

S(cid,speedid).t = t;
S(cid,speedid).F = F(FtoSave);

S(cid,speedid).COI = COI;

% Figure d: Coherence for 1 trial x
% Figure d: Coherence for 1 trial y
[WCOHx,WCSx,F,COI] = fCoherenceAnalysis(TRAJ(1,1).dx, TRAJ(1,1).dx, 9, 32);

% TBD: load trial 1

t = 1;
range = 100;
Freq = F(end-range:end);

c=1;
s=1;
subject = 1;

figure;
% p = panel();
% p.pack(2, 2); % subpanel
% %p.de.margin = 5;
% p.fontsize = 10;
% p.fontname = 'Arial';
% 
% % Org panels

subplot(1,2,1)
imagesc(TMP(c,s).COHx(:,:,n),[0 1]); colorbar;
title(sprintf('coh X: c %d s %d, F %2.2f\n', c,s, S(c,s).F));
hold on; plot(range-(S(c,s).COI.*(range/Freq(1))),'--w');
        
d = abs(Freq-S(c,s).F);
FtoSave = find(d==min(d));

set(gca, 'YTick',1:10:range,'YTickLabel',fix(Freq(1:10:range)*10)/10);
hold on; plot([0 4000],[FtoSave FtoSave],'--r','LineWidth',1);

subplot(1,2,2)

imagesc(TMP(c,s).COHy(:,:,n),[0 1]); colorbar;
title(sprintf('coh X: c %d s %d, F %2.2f\n', c,s, S(c,s).F));
hold on; plot(range-(S(c,s).COI.*(range/Freq(1))),'--w');
        
d = abs(Freq-S(c,s).F);
FtoSave = find(d==min(d));

set(gca, 'YTick',1:10:range,'YTickLabel',fix(Freq(1:10:range)*10)/10);
hold on; plot([0 4000],[FtoSave FtoSave],'--r','LineWidth',1);

figure
subplot(2,1,1)
plot(S(1,1).ointdx(n,:)'); hold on; plot(S(1,1).ttdx(n,:)');

subplot(2,1,2)
plot(S(1,1).ointdy(n,:)'); hold on; plot(S(1,1).ttdy(n,:),'r');
hold on; plot(S(1,1).ointdx(n,:)'); hold on; plot(S(1,1).ttdx(n,:),'r');


%%
figure; 
t = 1;
for c = 1:3
    for s = 1:3
        subplot(3,3,t);
        imagesc(circ_mean(G(c,s).GPHASEy,[],3),[-pi pi]); colorbar;
        title(sprintf('Phase Y: c %d s %d, F %2.2f\n', c,s, S(c,s).F));
        hold on; plot(range-(S(c,s).COI.*(range/Freq(1))),'--w');
        
        % target principal freq:         
        d = abs(Freq-S(c,s).F);
        FtoSave = find(d==min(d));

        set(gca, 'YTick',1:10:range,'YTickLabel',fix(Freq(1:10:range)*10)/10);
        hold on; plot([0 4000],[FtoSave FtoSave],'--r','LineWidth',1);

        t = t+1;
    end
end

figure;
%DSFigDef(2,[1. 1. 1],1);
t = 1;
for c = 1:3
    for s = 1:3
        subplot(3,3,t);
        imagesc(nanmean(G(c,s).GCOHx,3),[0 1]); colorbar;

        title(sprintf('COH X: c %d s %d, F %2.2f\n', c,s, S(c,s).F));
        hold on; plot(range-(S(c,s).COI.*(range/Freq(1))),'--w');
        
        % target principal freq:         
        d = abs(Freq-S(c,s).F);
        FtoSave = find(d==min(d));

        set(gca, 'YTick',1:10:range,'YTickLabel',fix(Freq(1:10:range)*10)/10);
        hold on; plot([0 4000],[FtoSave FtoSave],'--r','LineWidth',1);

        t = t+1;
    end
end


%% coherence and matches correlation
% Get number for coherence
t = 1;
for c =1:3
    for s = 1:3
        x = nanmean(G(c,s).COHx, 1);
        subx = nanmean(x,2);
        y = nanmean(G(c,s).COHy, 1);
        suby = nanmean(y,2);

        for n = 1:length(sub_list)
            Cx(t,:)=[subx(:,:,n) c s sub_list(n)];
            Cy(t,:)=[suby(:,:,n) c s sub_list(n)];

            t=t+1;
        end
    end
end
[sorted sid] = sort(Cx(:,4));
Cx = Cx(sid,:);
[sorted sid] = sort(Cy(:,4));
Cy = Cy(sid,:);

figure;
t = 1;
for c = 1:3
    for s = 1:3
        subplot(3,3,t)
        id = Cx(:,2)==c & Cx(:,3)==s;
        title(sprintf('c%d, s%d\n', c, s));
        plot(Cx(id,1),M(id,1),'o')
        corr([Cx(id,1),M(id,1)])
        t = t + 1;
    end
end
figure;
t = 1;
for c = 1:3
    for s = 1:3
        subplot(3,3,t)
        id = Cy(:,2)==c & Cx(:,3)==s;
        title(sprintf('c%d, s%d\n', c, s));
        plot(Cy(id,1),M(id,1),'o')
        corr([Cy(id,1),M(id,1)]);
        t = t + 1;
    end
end

%% Circular ANOVA ... Effect of cond, and vel.

%% HELP WITH INTERPRETATION: RUN THE WHOLE DEFAULT WCOHERENCE WITH A PURSUIT + TARGET + NOISE, COMPARING COIs ...
%%
% Pursuit signal: Simple ramp + steady, filtered (remove high-frequ comp) +
% white noise (wide spectrum) + target oscillation ...

% tx and ty
% all conditions
%len = length(TRAJ(1,1).dx);
%vel = mean(TRAJ(1,1).dx);

%latency =zeros(1,150);
%ramp = (0:(200-1)).*vel./(200-1);
%steady_state = ones(1,len-(length(latency)+length(ramp))).*vel;
%pursuit_ramp = [latency, ramp, steady_state]+TRAJ(c,s).dx'-vel+(rand(1,len).*4-2);

% a ramp filtered in time (note saccades are filtered ...; we can also simulate the effect of the interpolation!!)
        len = length(TRAJ(c,s).dx);
        vel = mean(TRAJ(c,s).dx);

        latency =zeros(1,150);
        ramp = (0:(200-1)).*vel./(200-1);   
        steady_state = ones(1,len-(length(latency)+length(ramp))).*vel;

pursuit_ramp = [latency, ramp, steady_state]+TRAJ(c,s).dx'-vel + (rand(1,len).*4-2);
figure; plot(pursuit_ramp)

        len = length(TRAJ(c,s).dy);
        vel = mean(TRAJ(c,s).dy);

        latency =zeros(1,150);
        ramp = (0:(200-1)).*vel./(200-1);   
        steady_state = ones(1,len-(length(latency)+length(ramp))).*vel;

pursuit_ramp = [latency, ramp, steady_state]+TRAJ(c,s).dy'-vel + (rand(1,len).*4-2);
figure; plot(pursuit_ramp)


figure
% x
t = 1;
for c = 1:3
    for s = 1:3
        subplot(3,3,t);
        
        len = length(TRAJ(c,s).dx);
        vel = mean(TRAJ(c,s).dx);

        latency =zeros(1,150);
        ramp = (0:(200-1)).*vel./(200-1);   
        steady_state = ones(1,len-(length(latency)+length(ramp))).*vel;
        pursuit_ramp = [latency, ramp, steady_state]+[latency TRAJ(c,s).dx(length(latency)+1:end)'-vel]+ (rand(1,len).*4-2);
        
        % Indetermination limits much larger!
        [wcoh,~,F,coi]=wcoherence(pursuit_ramp,TRAJ(c,s).dx',1000, 'NumOctaves', 10, 'VoicesPerOctave',32,'PhaseDisplayThreshold',.5);
        wcoherence(pursuit_ramp,TRAJ(c,s).dx',1000, 'NumOctaves', 10, 'VoicesPerOctave',32,'PhaseDisplayThreshold',1);
        %helperPlotCoherence(WCOH,tm,F,coi,'Seconds','Hz');
        
        t = t + 1;
        
    end
end

%%Y
figure
% x
t = 1;
for c = 1:3
    for s = 1:3
        subplot(3,3,t);
        
        len = length(TRAJ(c,s).dy);
        vel = mean(TRAJ(c,s).dy);

        latency =zeros(1,150);
        ramp = (0:(200-1)).*vel./(200-1);   
        steady_state = ones(1,len-(length(latency)+length(ramp))).*vel;
        pursuit_ramp = [latency, ramp, steady_state]+[latency TRAJ(c,s).dy(length(latency)+1:end)'-vel] + (rand(1,len).*6-3);
        
        % Indetermination limits much larger!
        [wcoh,~,F,coi]=wcoherence(pursuit_ramp,TRAJ(c,s).dy',1000, 'NumOctaves', 10, 'VoicesPerOctave',32,'PhaseDisplayThreshold',.5);
        wcoherence(pursuit_ramp,TRAJ(c,s).dy',1000, 'NumOctaves', 10, 'VoicesPerOctave',32,'PhaseDisplayThreshold',1);
        %helperPlotCoherence(WCOH,tm,F,coi,'Seconds','Hz');
        
        t = t + 1;
        
    end
end

