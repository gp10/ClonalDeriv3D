%% LOAD DATA ON CLONE PROJECTED AREAS:
load ClonalData_invitro.mat
mydirectory = pwd;

%% CALCULATE EXPERIMENTAL AVG. CLONE SIZE:
[mean_clonesizes,sem_clonesizes,mean_indiv_clonesizes,mean_All_clonesizes,std_All_clonesizes] = meanSizeCalc(rtime,nmice,rs,rs_all);

%% PLOT EXPERIMENTAL AVG. CLONE SIZE:
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
        plot(rtime,mean_indiv_clonesizes(ata,:),symbtype(ata),'Color',colWT(ata,:),'MarkerFaceColor',colWT(ata,:),'MarkerSize',5)
        % Fitting per mouse:
        LineFit = fit(rtime',mean_indiv_clonesizes(ata,:)','poly1');
        yfit = LineFit.p1*[0:30]+LineFit.p2;
        hold on
        plot([0:30],yfit,'--','Color',colWT(ata,:))
    end
end

% Display summary data:
for ete = 1:length(rtime)
    h = errorbar(rtime(ete),mean_clonesizes(1,ete),sem_clonesizes(1,ete),'CapSize',0); % SEM around avg clone size from all mice
    set(h,'Marker','.','MarkerSize',12,'Color','k','LineWidth',0.25)
end

xlim([3.5 20]); set(gca,'XTick',[4:2:20]);
ylim([0 2.5E5])
tobs = rtime;
xlabel('Time (weeks)')
ylabel('Average clone area (pixel^2)')
axis square

%% CALCULATE AND DISPLAY BEST FIT ON ALL AVG. CLONE SIZE DATA (from 3rd time point onwards):
myfun = fittype('a0+b*t','dependent',{'y'},'independent',{'t'},'coefficients',{'a0', 'b'});
LineFit = fit(reshape(repmat(rtime(3:end),nmice,1),nmice*(length(rtime)-2),1) - rtime(3), ...
    reshape(mean_indiv_clonesizes(:,3:end), size(mean_indiv_clonesizes(:,3:end),1)*size(mean_indiv_clonesizes(:,3:end),2),1), myfun);
LineFit_ci = confint(LineFit,0.95);
yfit = LineFit.b*[-1:30]+LineFit.a0;
% best fit on AvgCloneSize:
paramFit_avg = LineFit
save paramFit_avg.mat paramFit_avg
plot(([-1:30]+rtime(3)),yfit,'k--')

%% DISPLAY OPTIMUM PARAMETER FIT (obtained from fitting both Avg. Clone Size and Clone Density)
if exist('paramFit_optim.mat')
    load paramFit_optim.mat
    yfit2 = paramFit_optim.b*[-1:30]+paramFit_optim.a0;
    hold on
    plot(([-1:30]+rtime(3)),yfit2,'k-')
end

%% FUNCTIONS
% Average clone size calculation:
function [mean_size,sem_size,mean_size_indiv,mean_size_All,std_size_All] = meanSizeCalc(rtime,nmice,rx,rx_all)
    mean_size_indiv = [];
    mean_size = [];
    sem_size = [];
    for ata = 1:size(rtime,2)
        for ete = 1:nmice
            row_indiv_persis = find(rx{ete,1}(:,ata) > 0);
            mean_size_indiv(ete,ata) = mean(rx{ete,1}(row_indiv_persis,ata));
        end
        % Statistics in batch:
        row_All_persis = find(rx_all(:,ata) > 0);
        mean_size_All(1,ata) = mean(rx_all(row_All_persis,ata));
        std_size_All(1,ata) = std(rx_all(row_All_persis,ata),0,1);
    end
    % Statistics separating per animal:
    mean_size = mean(mean_size_indiv,1);
    sem_size = std(mean_size_indiv,0,1) ./ sqrt(nmice);
end

