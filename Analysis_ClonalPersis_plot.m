%% LOAD DATA ON CLONE PROJECTED AREAS:
load ClonalData_invitro.mat
mydirectory = pwd;

%% PLOT EXPERIMENTAL CLONAL DENSITY OVER TIME:
DisplayIndiv = false; % boolean (display data from individual experiments separately)
figure()
hold on

% Display individual-experiment data:
if DisplayIndiv
    % colorkey:
    colWT = [1 0 0;
        0.93 0.69 0.13;
        0.39 0.83 0.07;
        0.07 0.62 1;
        0.72 0.27 1;
        0 0 0];
    % symboltype:
    symbtype = 'osdvp*';
    for ata = 1:nmice
        plot(rtime, NClones(ata,:)./NClones(ata,1).*100 ,symbtype(ata),'Color',colWT(ata,:),'MarkerFaceColor',colWT(ata,:),'MarkerSize',5)
        %Fits per mouse: "No" will represent the initial population (unknown)
        myfun = fittype('No/(1+rlambda*t)','dependent',{'y'},'independent',{'t'},'coefficients',{'No', 'rlambda'});
        myfun_ini = [100 0.4]; myfun_lower = [5 0]; myfun_upper = [1000 Inf];
        [myfun_fit,myfun_fit_stat] = fit(rtime(2:end)' - rtime(2),NClones(ata,2:end)',myfun,'Startpoint',myfun_ini,'Lower',myfun_lower,'Upper',myfun_upper);
        infer_y = myfun_fit.No./(1+myfun_fit.rlambda.*[0:0.1:53]);
        plot([0:0.1:53] + rtime(2), infer_y./NClones(ata,1).*100,'--','Color',colWT(ata,:))
    end
end

% Calculate and display summary data:
mean_DensClones = mean(NClones./NClones(:,1));
sem_DensClones = std(NClones./NClones(:,1),0,1) ./ sqrt(nmice);
h1 = errorbar(rtime,mean_DensClones.*100,sem_DensClones.*100,'CapSize',0); % around avg survival per mice
set(h1,'LineStyle','none','Marker','.','MarkerSize',12,'Color','k','LineWidth',0.25)

xlim([3.5 20]); set(gca,'XTick',[4:2:20]);
ylim([0 105])
%tobs = rtime;
xlabel('Time (weeks)')
ylabel('Number of clones (%)')
axis square

%% CALCULATE AND DISPLAY BEST FIT ON ALL CLONAL DENSITY DATA (from 3rd time point onwards):
myfun = fittype('No/(1+decay*t)','dependent',{'y'},'independent',{'t'},'coefficients',{'No', 'decay'});
myfun_ini = [100 0.4]; myfun_lower = [0 0]; myfun_upper = [1000 Inf];
[myfun_fit,myfun_fit_stat] = fit(reshape(repmat(rtime(3:end),nmice,1),nmice*(length(rtime)-2),1) - rtime(3), ...
    reshape(NClones(:,3:end)./NClones(:,1), size(NClones,1)*size(NClones(:,3:end),2), 1), ...
    myfun,'Startpoint',myfun_ini,'Lower',myfun_lower,'Upper',myfun_upper);
infer_y = myfun_fit.No./(1+myfun_fit.decay.*[-0.5:1:53]);
% best fit on Clonal Density decay:
paramFit_decay = myfun_fit
save paramFit_decay.mat paramFit_decay
plot(([-0.5:1:53] + rtime(3)), infer_y.*100,'k--')

%% DISPLAY OPTIMUM PARAMETER FIT (obtained from fitting both Avg. Clone Size and Clone Density)
if exist('paramFit_optim.mat')
    load paramFit_optim.mat
    infer_y2 = paramFit_optim.N0./(1+(paramFit_optim.b/paramFit_optim.a0).*[-0.5:1:53]);
    hold on
    plot(([-0.5:1:53] + rtime(3)), infer_y2.*100,'k-')
end

%% OTHER VALUES
% decay rate (b/a0) given by best-fit to avg-clone growth (a0+bt), where a0 is initial avg clone area and b is avg clone growth rate, so that total population remains constant as a0*No.
if exist('paramFit_avg.mat')
    load paramFit_avg.mat
    decay_inferred = (paramFit_avg.b/paramFit_avg.a0)
end
