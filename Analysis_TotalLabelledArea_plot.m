%% LOAD DATA ON CLONE PROJECTED AREAS:
load ClonalData_invitro.mat
mydirectory = pwd;

%% COMPUTE TOTAL LABELLED AREA PER EXPERIMENT:
AreaLabelled = [];
for aja = 1:nmice
    for eje = 1:length(rtime)
        AreaLabelled(aja,eje) = sum(rs{aja,1}(find(rs{aja,1}(:,eje)>0),eje),1);
    end
end

%% PLOT EXPERIMENTAL TOTAL LABELLED AREA OVER TIME:
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
        plot(rtime, AreaLabelled(ata,:),symbtype(ata),'Color',colWT(ata,:),'MarkerFaceColor',colWT(ata,:),'MarkerSize',5)
    end
end

% Calculate and display summary data:
mean_AreaLabelled = mean(AreaLabelled,1);
sem_AreaLabelled = std(AreaLabelled,0,1) ./ sqrt(nmice);
h1 = errorbar(rtime,mean_AreaLabelled, sem_AreaLabelled,'CapSize',0);  set(h1,'LineStyle','none','Marker','.','MarkerSize',12,'Color','k','LineWidth',0.25)

xlim([3.5 20]); set(gca,'XTick',[4:2:20]);
xlabel('Time (weeks)')
ylabel('Total labelled area'); ylim([0 4.5E6])
axis square

%% CALCULATE AND DISPLAY BEST FIT ON ALL TOTAL LABELLED AREA (from 3rd time point onwards):
myfun = fittype('a*t+b','dependent',{'y'},'independent',{'t'},'coefficients',{'a', 'b'});
myfun_ini = [0 1E4]; myfun_lower = [0 0]; myfun_upper = [0 Inf];
[myfun_fit,myfun_fit_stat] = fit(reshape(repmat(rtime(3:end),nmice,1),nmice*(length(rtime)-2),1) - rtime(3), ...
    reshape(AreaLabelled(:,3:end), size(AreaLabelled,1)*size(AreaLabelled(:,3:end),2), 1), ...
    myfun,'Startpoint',myfun_ini,'Lower',myfun_lower,'Upper',myfun_upper);
infer_y = myfun_fit.a.*[-1:1:53] + myfun_fit.b;
% best fit on Total Labelled Area:
plot(([-1:1:53] + rtime(3)), infer_y,'k--')

%% DISPLAY OPTIMUM PARAMETER FIT (obtained from fitting both Avg. Clone Size and Clone Density):
if exist('paramFit_optim.mat')
    load paramFit_optim.mat
    % convert N0 fraction into absolute value:
    infer_y2 = zeros.*[-1:1:53] + paramFit_optim.a0 * ( paramFit_optim.N0 * mean(NClones(:,1)) );
    hold on
    plot(([-1:1:53] + rtime(3)), infer_y2,'k-')
end
