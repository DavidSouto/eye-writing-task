%%%%%%%%%%%%%%%%%%%% (EL2_playback: courtesy of Alexander C Schutz @ Giessen) %%%%%%%%%%%
%['\data\e',int2str(trl.exp),'s',int2str(trl.sub),
%util.el_datname]
status = EL2_playback(util.el_datdir);
if ~status;
    error('playback failed');
else

    sLoadOnline;                        

end                
                    