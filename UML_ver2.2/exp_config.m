function par = exp_config()

% set the experimental parameters, the parameter space and the priors to 
% be used as an input argument of the UML class. 
% For example, uml = UML(par).

par.model = 'logit';    % the model for the psychometric function, choose 
                        % between 'logit' for the logistic psychometric
                        % function, or 'weibull' for the Weibull
                        % psychometric function.
                        
% The logistic psychometric function is given by:
%
% p = gamma+(1-gamma-lambda).*(1+exp(-(x-alpha).*beta)).^(-1),
%
% where p is proportion correct, x is the signal strength, gamma indicates 
% chance performance (e.g., gamma = 0.5 for 2AFC experiments), lambda is 
% the lapse rate (difference between the upper asymptote of the 
% psychometric function and 1, alpha and beta are the threshold and beta 
% parameters, respectively.
%
% The Weibull psychometric fucntion is given by:
%
% p = gamma+(1-gamma-lambda).*(1-exp(-(x./alpha).^beta)).
%
% Note that, for Weibull psychometric function, the alpha parameter has to
% be strictly greater than zero, and the beta parameter has to be greater 
% or equal to one.

par.ndown = 2;  % the parameter for the up-down sweetpoint selection rule
par.method = 'mean';    % the method for estimating the parameters from the
                        % the posterior parameter distribution. choose
                        % between 'mean' and 'mode'.
par.x0 = 30;    % the initial signal strength
par.x_lim = [-30 30];   % the limits to the signal strength

par.alpha = struct(...
    'limits',[-10 30],...       %range of the parameter space for alpha
    'N',61,...                %number of alpha values. If this value is set to 1, then the first element of alpha.limits would be the assumed alpha and the alpha parameter is not estimated.
    'scale','lin',...         %the linear or log spacing. Choose between 'lin' and 'log'.
    'dist','flat',...         %prior distribution of the alpha parameter. Choose between 'norm' and 'flat'.
    'mu',0,...                %mean of the prior distribution.
    'std',20 ...              %standard deviation of the prior distribution.  
    );

par.beta = struct(...
    'limits',[.1 10],...      %range of the parameter space for beta
    'N',41,...                %number of beta values. If this value is set to 1, then the first element of beta.limits would be the assumed beta and the beta parameter is not estimated.
    'scale','log',...         %the linear or log spacing. Choose between 'lin' and 'log'.
    'dist','flat',...         %prior distribution of the beta parameter. Choose between 'norm' and 'flat'.
    'mu',1,...                %mean of the prior distribution.
    'std',5 ...               %standard deviation of the prior distribution.
    );

par.gamma = 0.5;

par.lambda = struct(...
    'limits',[0 0.2],...      %range of the parameter space for lambda
    'N',11,...                 %number of lambda values. If this value is set to 1, then the first element of lambda.limits would be the assumed lambda and the lambda parameter is not estimated.
    'scale','lin',...         %the linear or log spacing. Choose between 'lin' and 'log'.
    'dist','flat',...         %prior distribution of the lambda parameter. Choose between 'norm' and 'flat'.
    'mu',0,...                %mean of the prior distribution.
    'std',0.1 ...             %standard deviation of the prior distribution.  
    );

% Note that, when 'scale' is set to 'lin' and 'dist' is set to 'norm', the 
% prior distribution is constructed based on:
% p0 ~ N(mu,std)
% On the other hand, when 'scale' is set to 'lin' and 'dist' is set to 
% 'norm', the prior distribution is constructed based on:
% log10(p0) ~ N(mu,std)

end

%eof