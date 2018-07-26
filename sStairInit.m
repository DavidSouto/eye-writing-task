% stair input
% stair input
sOpt.level_set    = -10:0.2:10;       % set of possible probe values
sOpt.pse_set     = -10:0.2:10;       % sampling of pses, doesn't have to be the same as probe set
sOpt.slope_set    = (.5:.1:5).^2;     % set of slopes, quad scale
sOpt.lapse       = 0.05;             % lapse / mistake rate
sOpt.guess       = 0;             % guess rate / minimum correct response rate (for detection expt: 1/num_alternative, 0 for discrimination expt)

% general settings
ntrial  = trl.STAIRCASE_NTRIALS;

stair = MinExpEntStair('v2');
% use lookup table to cache pvalues and avoid unnecessary evaluations of
% psychometric function? Can require lots of memory, especially when
% stepsize of probeset and meanset is not equal. Call before calling
% stair.init.
stair.set_use_lookup_table(true);
% option: use logistic instead of default cumulative normal. best to call
% before stair.init
% stair('set_psychometric_func','logistic');
% init stair
stair.init(sOpt.level_set,sOpt.pse_set,sOpt.slope_set,sOpt.lapse,sOpt.guess);
% option: use a subset of all data for choosing the next probe, use
% proportion of available data (good idea for robustness - see docs)
stair.toggle_use_resp_subset_prop(10,.9);

% option: instead of choosing a random value for the first probe,
% you can set which value is to be tested first.
stair.set_first_value(0);

