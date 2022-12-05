%% LOAD DATA ON CLONE PROJECTED AREAS:
load ClonalData_invitro.mat
mydirectory = pwd;

% Include avg clone size calculations:
[mean_clonesizes,sem_clonesizes,mean_indiv_clonesizes,mean_All_clonesizes,std_All_clonesizes] = meanSizeCalc(rtime,nmice,rs,rs_all);

%% CALCULATE GoF VALUES IN A PARAMETER GRID SEARCH:
% Parameters to search for:
all_b = [6e3:1e1:1.6e4];
all_a0 = [3e4:1e2:4.5e4];
all_N0 = [0.72:0.01:0.90];

% Iteration on GoF calculation
GoF_avg = []; GoF_decay = []; GoF_overall = [];
for iter = 1:length(all_b)
    b = all_b(iter);
    for itir = 1:length(all_a0)
        a0 = all_a0(itir);
        for itor = 1:length(all_N0)
            N0 = all_N0(itor);

            % GoF on average clone size:
            yfit = b*(rtime(3:end)-rtime(3)) + a0; % avgCloneSize should scale linearly with time
            subset_obs = mean_indiv_clonesizes(:,3:end);
            GoF_avg(iter,itir,itor) = GoF(subset_obs,yfit);

            % GoF on clone density:
            yfit = N0 ./(1+ (b/a0).*(rtime(3:end)-rtime(3)) ); % clone density should scale with time following an hyperbolic decay
            subset_obs = NClones(:,3:end)./NClones(:,1);
            GoF_decay(iter,itir,itor) = GoF(subset_obs,yfit);

            % Overall GoF:
            GoF_overall(iter,itir,itor) = mean( [GoF_avg(iter,itir,itor) GoF_decay(iter,itir,itor)] );
        end
    end
end

% Optimum parameter values (based on GoF):
[xloc,yloc,zloc] = ind2sub([length(all_b) length(all_a0) length(all_N0)], find(GoF_overall == min(GoF_overall(:))) );
paramFit_optim.b = all_b(xloc);
paramFit_optim.a0 = all_a0(yloc);
paramFit_optim.N0 = all_N0(zloc);
paramFit_optim
save paramFit_optim.mat paramFit_optim

%% PLOT GoF VALUES IN THE PARAMETER SPACE:
figure()
for itor = 1:length(all_N0)
    ax = subplot(4,5,itor);
    subdata(:,:) = GoF_overall(:,:,itor);
    imagesc(subdata,[min(GoF_overall(:)) max(GoF_overall(:))])
    ylabel('b')
    set(gca,'YTick',[1 length(all_b)])
    set(gca,'YTickLabel',num2cell([all_b(1), all_b(end)]))
    set(gca,'YDir','normal')
    xlabel('a0')
    set(gca,'XTick',[1 length(all_a0)])
    set(gca,'XTickLabel',{sprintf('%.1e',all_a0(1)), sprintf('%.1e',all_a0(end))})
    title(sprintf(['N_0 = %.2f'], all_N0(itor)))
end
subplot(4,5,20); axis off
colorbar

% Plot optimum parameter values:
subplot(4,5,zloc)
hold on
plot(yloc,xloc,'r*')

%% FUNCTIONS
% Goodness-of-Fit (GoF) statistic calculation:
function GoF_val = GoF(yobs,yfit)
    mean_obs = mean(yobs(:));
    Stot = sum( (yobs(:) - mean_obs).^2 );
    Sres = 0;
    for myt = 1:size(yobs,2)
        Sres = Sres + sum( (yobs(:,myt) - yfit(1,myt)).^2 );
    end
    GoF_val = Sres / Stot;
end

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
