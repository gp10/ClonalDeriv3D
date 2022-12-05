%% RUN -OR- LOAD SIMULATION DYNAMICS:
% RUN simulation if and only if simulation outcome has not been stored yet:
if ~exist('./Simulation_data/MCMC_a_1000_Dim_200_t_26w_optim.mat')
    
    %% ESTIMATE OF GRID-CELL REPLACEMENT RATE (Lambda):
    % Set realistic values:
    if exist('paramFit_optim.mat')
        load paramFit_optim.mat
    else
        Analysis_OptimumFit
    end
    
    % Set Area corresponding to grid unitary cell:
    %a_gridcell = min(rs_all(:)); % pixel^2 / gridcell
    a_gridcell = 1000; % pixel^2 / gridcell
    
    % Area growth rate:
    %paramFit_avg.b; % pixel^2 / time
    paramFit_optim.b; % pixel^2 / time
    
    % Gridcell-based growth rate:
    %b_gridcell = paramFit_avg.b / a_gridcell % gridcell / time
    b_gridcell = paramFit_optim.b / a_gridcell % gridcell / time
    
    % Gridcell-based replacement rate:
    Lambda = b_gridcell/0.5; % each progenitor cell can replace or be replaced with equal chance (50%) per replacement event.
    
    %% OTHER PARAMETER SETTINGS FOR SIMULATION RUN:
    timelim = 26; % Time interval
    freqLabel = 0.15; % Frequency of initially labelled cells
    lattice.Dim = 200; % Lattice size (width)
    lattice.Neigh = 6; % No. of neighbors
    nval = 260; % No. of timepoints where to retrieve data from
    indiv = 1; % No. of independent simulations

    %% RUN SIMULATION OF IN VITRO CLONE DYNAMICS IN A 2D GRID:
    [nx_basal,ntime,ALL_x_Clone,ALL_x_Label] = Simul_2Dgrid_SPdynamics(timelim,Lambda,freqLabel,lattice,nval,indiv);
    
else
    %% LOAD SIMULATION RESULTS:
    load ./Simulation_data/MCMC_a_1000_Dim_200_t_26w_optim.mat
end

%% EXTRACT AND SHOW STATISTICAL PROPERTIES OF SIMULATED CLONE SIZES:
isizes(:,:) = nx_basal(1,:,:);

% Check surviving clone grid sizes at end timepoint:
isizes(end, find(isizes(end,:)>0) )

% Store inferred average clone grid sizes at different timepoints:
misavg_gridcell = [];
for aja = 1:nval+1
    misavg_gridcell(aja) = mean ( isizes(aja, find(isizes(aja,:)>0) ) );
end

% Convert average clone grid sizes into average clone areas:
misavg_a = misavg_gridcell * a_gridcell;

% Check inferred average clone area at end timepoint:
mean ( isizes(end, find(isizes(end,:)>0) ) ) * a_gridcell

% Store inferred clone density at different timepoints:
misdens = [];
for aja = 1:nval+1
    misdens(aja) = length(find(isizes(aja,:)>0));
end

% Set time-delay vs. experimental data so as to fit AvgCloneSize at timepoint 3 (initial condition):
load ClonalData_invitro.mat
ref_t0_pos = find(misavg_a >= paramFit_optim.a0,1);
tdelay = rtime(3) - ntime( ref_t0_pos ) % weeks

% Scale simulated No. surviving clones to fit relative clone density at timepoint 3 (initial condition):
scale_factor = misdens(ref_t0_pos) / paramFit_optim.N0; % absolute-to-relative value conversion
misdens_scaled = misdens / scale_factor;
%paramFit_optim.a0 / scale_factor

%% GENERATE SUBSAMPLES OF SIMULATED CLONES GIVEN EXPERIMENTAL SAMPLE SIZE:
% Recast statistical properties of experimental sample size:
mean_DensClones = mean(NClones./NClones(:,1));
sem_DensClones = std(NClones./NClones(:,1),0,1) ./ sqrt(nmice);

% Experimental sample size at 3rd time point (taken as initial time point for simulations):
t0_NClones_mean = mean_DensClones(3) * mean(NClones(:,1));
t0_NClones_sem = sem_DensClones(3) * mean(NClones(:,1));

% Simulated clone size subsampling
isizes_samples = {};
n_samples = 500; % No. of subsamples
for iter = 1:n_samples
    t0_NClones_sample = round( randn(1)*t0_NClones_sem + t0_NClones_mean ); % random pick from a normal dist with SD equal to experimental sem (convert from z-dist)
    loc_t0_survClones = find(isizes(ref_t0_pos,:)>0);
    isizes_samples{iter,1} = isizes(:,loc_t0_survClones( randperm( length(loc_t0_survClones), t0_NClones_sample ) ));
end

%% FIT SIMULATED AVERAGE CLONE SIZE TO EXPERIMENTAL AVERAGE CLONE SIZE:
% Experimental average clone size:
Analysis_AvgCloneSize_plot

% Compute Avg clone size in simulated subsamples:
sim_avgsize = [];
for iter = 1:n_samples
    for itir = 1:nval+1
        sim_avgsize(iter,itir) = mean ( isizes_samples{iter,1}(itir, find(isizes_samples{iter,1}(itir,:)>0)) );
    end
end
sim_avgsize_a = sim_avgsize * a_gridcell;

sim_avgsize_a_mean = mean(sim_avgsize_a,1);
sim_avgsize_a_std = std(sim_avgsize_a,0,1);

% Plot simulated average clone size:
hold on
plot((ntime+tdelay),misavg_a, 'Color', [1 0.41 0.16]) % overall
%plot((ntime(ref_t0_pos:end)+tdelay),sim_avgsize_a_mean(ref_t0_pos:end), 'Color', [1 0.41 0.16])

% Plot simulated plausible interval on average clone size:
hold on
plot((ntime(ref_t0_pos:end)+tdelay),sim_avgsize_a_mean(ref_t0_pos:end)+sim_avgsize_a_std(ref_t0_pos:end), 'Color', [1 0.41 0.16])
plot((ntime(ref_t0_pos:end)+tdelay),sim_avgsize_a_mean(ref_t0_pos:end)-sim_avgsize_a_std(ref_t0_pos:end), 'Color', [1 0.41 0.16])

%% FIT SIMULATED CLONE DENSITY TO EXPERIMENTAL CLONE DENSITY:
% Experimental clone density:
Analysis_ClonalPersis_plot

% Compute Clone density in simulated subsamples:
sim_NClones = [];
for iter = 1:n_samples
    for itir = 1:nval+1
        sim_NClones(iter,itir) = length( find(isizes_samples{iter,1}(itir,:) > 0) );
    end
end
sim_NClones_avg = mean(sim_NClones,1);
sim_NClones_up = mean(sim_NClones,1) + std(sim_NClones,0,1);
sim_NClones_dn = mean(sim_NClones,1) - std(sim_NClones,0,1);

% Plot simulated clone density:
hold on
plot((ntime+tdelay),(misdens_scaled)*100, 'Color', [1 0.41 0.16]) % take cellsize as smallest a_gridcell value (1,000 pixel^2)
%plot((ntime(ref_t0_pos:end)+tdelay),(sim_NClones_avg(ref_t0_pos:end)/t0_NClones_mean)*paramFit_optim.N0*100, 'Color', [1 0.41 0.16]) % take cellsize as smallest a_gridcell value (1,000 pixel^2)

% Plot simulated plausible interval on clone density:
hold on
plot((ntime(ref_t0_pos:end)+tdelay),(sim_NClones_up(ref_t0_pos:end)/t0_NClones_mean)*paramFit_optim.N0*100, 'Color', [1 0.41 0.16]) % take cellsize as smallest a_gridcell value (1,000 pixel^2)
plot((ntime(ref_t0_pos:end)+tdelay),(sim_NClones_dn(ref_t0_pos:end)/t0_NClones_mean)*paramFit_optim.N0*100, 'Color', [1 0.41 0.16]) % take cellsize as smallest a_gridcell value (1,000 pixel^2)
% fill([(ntime(ref_t0_pos:end)+tdelay) fliplr(ntime(ref_t0_pos:end)+tdelay)], ...
%     [(sim_NClones_up(ref_t0_pos:end)/t0_NClones_mean) fliplr(sim_NClones_dn(ref_t0_pos:end)/t0_NClones_mean)]*paramFit_optim.N0*100, ...
%     [1 0.41 0.16]); % take cellsize as smallest a_gridcell value (1,000 pixel^2)

%% FIT SIMULATED TOTAL LABELLED AREA TO EXPERIMENTAL TOTAL LABELLED AREA:
% Experimental total labelled area:
Analysis_TotalLabelledArea_plot

% Compute total labelled area in simulated subsamples:
sim_totallabeled = [];
for iter = 1:n_samples
    for itir = 1:nval+1
        sim_totallabeled(iter,itir) = sum ( isizes_samples{iter,1}(itir, find(isizes_samples{iter,1}(itir,:)>0)) );
    end
end
sim_totallabeled_a = sim_totallabeled * a_gridcell;

sim_totallabeled_a_mean = mean(sim_totallabeled_a,1);
sim_totallabeled_a_std = std(sim_totallabeled_a,0,1);

% Plot simulated total labelled area:
hold on
plot((ntime(ref_t0_pos:end)+tdelay),sim_totallabeled_a_mean(ref_t0_pos:end), 'Color', [1 0.41 0.16])

% Plot simulated plausible interval on total labelled area:
hold on
plot((ntime(ref_t0_pos:end)+tdelay),sim_totallabeled_a_mean(ref_t0_pos:end) + ...
    sim_totallabeled_a_std(ref_t0_pos:end), 'Color', [1 0.41 0.16])
plot((ntime(ref_t0_pos:end)+tdelay),sim_totallabeled_a_mean(ref_t0_pos:end) - ...
    sim_totallabeled_a_std(ref_t0_pos:end), 'Color', [1 0.41 0.16])

%%
% myvector = [1000 2000 5000 8000 10000 20000];
% tdelay = [42  42.44  43.78  45.12  46.01  50.73];
% mycol = 'rygbmk';
% %figure()
% for eje = 1:length(myvector)
%     eval(sprintf("load MCMC_a_%d_Dim_200_t_10w.mat", myvector(eje)))
%     Lambda*0.5 * a_gridcell
%     isizes(:,:) = nx_basal(1,:,:);
%     hold on
%     misavg = [];
%     for aja = 1:nval+1
%         misavg(aja) = mean ( isizes(aja, find(isizes(aja,:)>0) ) );
%     end
%     plot(ntime*7+tdelay(eje),misavg* a_gridcell, mycol(eje))
% end
% 
% figure()
% for eje = 1:length(myvector)
%     eval(sprintf("load MCMC_a_%d_Dim_200_t_10w.mat", myvector(eje)))
%     Lambda*0.5 * a_gridcell
%     isizes(:,:) = nx_basal(1,:,:);
%     isizes(end, find(isizes(end,:)>0) );
%     mean ( isizes(end, find(isizes(end,:)>0) ) ) * a_gridcell;
%     hold on
%     misdens = [];
%     for aja = 1:nval+1
%         misdens(aja) = length(find(isizes(aja,:)>0));
%     end
%     misdens
%     plot(ntime*7+tdelay(eje),(misdens / a_gridcell)*1000/40000*100, mycol(eje)) % take cellsize as smallest a_gridcell value (1,000 pixel^2)
% end
