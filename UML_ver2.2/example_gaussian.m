clc;
close all;
clear all;

uml = UML(exp_config_gaussian());
uml.setPhi0([-10,8,0.5,0.05]);

ntrials = 200;
for i = 1:ntrials 
    % present the stimulus and collect the observer's response r in in
    % terms of correct (1) or incorrect (0).
    r = uml.simulateResponse(uml.xnext);
    
    % update the signal level
    uml.update(r);
    
    uml.plotP();
end

figure;
subplot('position',[0.1 0.7 0.85 0.25],'fontsize',12);
plot(uml.phi(:,1),'k','LineWidth',2);hold;
line([0 ntrials],uml.phi0(1)*ones(1,2),'Color','k','LineStyle','--');
ylabel('\alpha');
axis([0 ntrials uml.par.alpha.limits+[-3 3]]);
set(gca,'XTickLabel',[]);

subplot('position',[0.1 0.4 0.85 0.25],'fontsize',12);
semilogy(uml.phi(:,2),'k','LineWidth',2);hold;
line([0 ntrials],uml.phi0(2)*ones(1,2),'Color','k','LineStyle','--');
ylabel('\beta');
axis([0 ntrials uml.par.beta.limits.*[0.5 2]]);
set(gca,'XTickLabel',[],'YTick',[0.1 1 10]);

subplot('position',[0.1 0.1 0.85 0.25],'fontsize',12);
plot(uml.phi(:,4),'k','LineWidth',2);hold;
line([0 ntrials],uml.phi0(4)*ones(1,2),'Color','k','LineStyle','--');
xlabel('Trial number');
ylabel('\lambda');
axis([0 ntrials uml.par.lambda.limits+[-0.05 0.05]]);