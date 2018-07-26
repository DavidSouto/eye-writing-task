[m,s,loglik]    = stair.get_fit();

[ps,rs]         = stair.get_history();

% Save stair structure

try
    %QuestBetaAnalysis(q); % optional
    save(['thres.' int2str(trl.sub) '.' int2str(trl.blk) '.mat'], 'stair','trl');

%     figure(1);
%     subplot(1,2,1);
%     imagesc(exp(.5*loglik));
%     set(gca,'YTick',1:4:length(sOpt.slope_set));
%     set(gca,'YTickLabel',sOpt.slope_set(1:4:end));
%     set(gca,'XTick',1:5:length(sOpt.pse_set));
%     set(gca,'XTickLabel',sOpt.pse_set(1:5:end));
%     title('estimated likelihood function');
%     xlabel('mean (a)')
%     ylabel('slope (b)')
% 
    subplot(1,2,2);
    pind = find(rs>0);
    nind = setdiff(1:length(ps),pind);
    plot(1:length(ps),ps,'k-',pind,ps(pind),'bo',nind,ps(nind),'ro');
    ylim([min(sOpt.level_set) max(sOpt.level_set)]);
    title('probe-resp history');

    [PSEfinal,DLfinal,loglikfinal]  = stair.get_PSE_DL();

    finalent = sum(-exp(loglikfinal(:)).*loglikfinal(:));
    fprintf('final estimates:\nPSE: %f\nDL: %f\nent: %f\n',PSEfinal,DLfinal,finalent);

    % for actual offline fitting of your data, you would probably want to use a
    % dedicated toolbox such as Prins, N & Kingdom, F. A. A. (2009) Palamedes:
    % Matlab routines for analyzing psychophysical data.
    % http://www.palamedestoolbox.org.
    % Also note that while the staircase runs far more rebust when a small
    % lapse rate is assumed, it is common to either fit the psychometric
    % function without a lapse rate, or otherwise with the lapse rate as a free
    % parameter (possibily varying only over subjects, but not over conditions
    % within each subject).

    switch stair.get_psychometric_func()
        case 'cumGauss'
            title('estimated likelihood function - cumulative Gaussian')
            ylabel('$\sigma$','interpreter','latex')
        case 'logistic'
            title('estimated likelihood function - logistic')
            ylabel('$s$','interpreter','latex')
    end
catch
    fprintf('Unable to save data\n');
    sca
    keyboard
    
end
