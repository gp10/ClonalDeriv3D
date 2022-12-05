%% RETRIEVE SIMULATION DATA
load ./Simulation_data/MCMC_a_1000_Dim_200_t_26w_optim.mat
Analysis_Simul_2Dgrid

%% COLOR SETTINGS
myColLabel = [0 0 0; 1 1 0]; % black or yellow
Xfrac = 0.3; % fraction of initial clones (cells) to label
myColConfetti = [0 1 0; ... % GFP
                 1 0 0; ... % RFP
                 1 1 0; ... % YFP
                 0 0 1];    % CFP
myColConfetti_frac = [60 20 5 1]; % proportion of each color relative to CFP
% Assign random colors to clones:
myColClones = myColConfetti( sum( mnrnd(1,myColConfetti_frac./sum(myColConfetti_frac),lattice.Dim^2) .* [1 2 3 4], 2), :);
myColClones = rand(lattice.Dim^2,3);
myColLabelClones = [0 0 0;...
            myColClones];
% Omit color for non-labelled clones:
myColClones = myColClones .* binornd(1,Xfrac,size(myColClones,1),1);

%% PLOTTING SLOTS OF 2D BASAL-LAYER SECTIONS
figure()
mytime = [0 4 8 12 16]; % selected simulation time to plot (weeks)
for aja = 1:length(mytime)
    mypoint = find(ntime >= mytime(aja),1);
    %x_Label_to_plot(:,:) = ALL_x_Label(1,mypoint,:,:);
    x_Clone_to_plot(:,:) = ALL_x_Clone(1,mypoint,:,:);
    %x_LabelClone_to_plot(:,:) = x_Label_to_plot .* x_Clone_to_plot;
    
    subplot(3,3,aja)
    imagesc(x_Clone_to_plot,[1 lattice.Dim^2])
    colormap(myColClones);
    freezeColors;
    ax = gca; axis(ax,'off'); axis square; 
    set(ax,'xaxisLocation','top')
    xlabel(ax,[sprintf('%0.1f',ntime(1,mypoint)+tdelay) ' weeks']); set(get(ax,'XLabel'),'Visible','on')
    %text(10,-10,['time = ' sprintf('%0.1f',ntime(1,mypoint)) ' weeks'],'Color','k','FontSize',12)
end

%% PLOTTING TIME EVOLUTION OF CELL FATES ON 2D BASAL-LAYER SECTION
figure('Color',[1 1 1])
hold on
for es = 1:length(ntime)
    cla
    %x_Label_to_plot(:,:) = ALL_x_Label(1,es,:,:);
    x_Clone_to_plot(:,:) = ALL_x_Clone(1,es,:,:);
    %x_LabelClone_to_plot(:,:) = x_Label_to_plot .* x_Clone_to_plot;
    
    %subplot(5,5,[1 25])
    imagesc(x_Clone_to_plot,[1 lattice.Dim^2])
    colormap(myColClones);
    freezeColors;
    ax = gca; axis(ax,'off'); axis square; 
    set(ax,'xaxisLocation','top')
    %xlabel(ax,['time = ' sprintf('%0.1f',ntime(1,es)) ' weeks']); set(get(ax,'XLabel'),'Visible','on')
    text(10,-10,['time = ' sprintf('%0.1f',ntime(1,es)+tdelay) ' weeks'],'Color','k','FontSize',12)
    axis off

    drawnow;
    mov(es)=getframe(gcf);
end
%         movie2avi(mov,'time_evol_basal_distrib.avi','compression','None');
vfile = VideoWriter('Video_2D_clones.avi');
%vfile.FileFormat = 'mp4';
vfile.Quality = 100; vfile.FrameRate = 24;
%vfile.VideoCompressionMethod = 'H.264'
%vfile.Height = 276; vfile.Width = 492;
open(vfile)
writeVideo(vfile,mov)
close(vfile)
