function[] = rungetgeometry_v2()
    % calls getgeometry_v2
    
    % set control parameters:
    t1 = now;                               % for timer
    diary 'getgeometry controlwindow output.txt'; diary on; % output the following to the diary
    filename = '.txt';  % filename with '.txt' suffix
    subsamplefreq = 20;                     % number of subsamples per recorded sample (total number of iterations are subsamplefreq*MCmax)
    nochains = 5;                           % number of independent chains
    sorttype = 'circumference';             % criterion when sorting state labels; 'none', 'circumference', 'proportion'
    [ mydata, nostates,~,~,nodata ] = readpolygonstates( filename );
    posstatevskin = true(nostates,nodata);  % unless declared otherwise, all measurements can be in all states
    %posstatevskin(2,(0+1):(100)) = false;  % first 100 measurements are informing dataset for first state
    datasetvskin = true(1,nodata); 
    %datasetvskin = false(2,nodata); datasetvskin(1,1:100) = true; datasetvskin(2,(100+1):end) = true;  % two subsets: First 100 measurements in one subset, all other ones in second subset
    
    fprintf( ' Info - rungetgeometry_v2: Start %d independent chain(s) now.\n', nochains );
    truegeo{nochains} = [];                 % declare
    currentclock = clock();
    parfor j_chain = 1:nochains
        truegeo{j_chain} = getgeometry_v2;
        truegeo{j_chain}.globaltimestamp = currentclock;
        truegeo{j_chain}.inputfilename = filename;
        truegeo{j_chain}.flpos_data = mydata;
        truegeo{j_chain}.MCmax = 1e4;                               % total number of recorded samples
        truegeo{j_chain}.burnin = ceil(0.4*truegeo{j_chain}.MCmax); % burnin (of recoreded samples)
        truegeo{j_chain}.subsamplefreq = subsamplefreq;
        truegeo{j_chain}.nostates = nostates;
        truegeo{j_chain}.posstatevskin = posstatevskin;
        truegeo{j_chain}.datasetvskin = datasetvskin;
        truegeo{j_chain}.nodsets = size(datasetvskin,1);
        truegeo{j_chain}.comment = sprintf('%d',j_chain);
        truegeo{j_chain}.initialiseparameters();
        truegeo{j_chain}.control();
    end     % end of chains loop
    fprintf( ' Info - rungetgeometry_v2: Final statistics (after %1.3f sec):\n', (now-t1)*3600*24 );
    getfinalstatistics( truegeo, sorttype );
    fprintf( ' Info - rungetgeometry_v2: Done sampling now (after %1.3f sec).\n', (now-t1)*3600*24 );
    diary 'getgeometry controlwindow output.txt'; diary off; % stop outputting to the diary
end     % end of rungetgeometry_v2 function

function[ mydata, nostates,nomarkers,dim,nodata ] = readpolygonstates( filename )
    
    % open file to read:
    skip = 0; maxnodata = Inf;                              % can be adjusted, if only subsets shall be analysed
    textfile = fopen( filename, 'r' );                      % open textfile to read
    fprintf( ' Info - readpolygonstates: Start reading:\n  %s\n', filename );
    fprintf( '  within [%d..%d]\n', 1,maxnodata );
    fgetl(textfile);    % date
    buffer = fgetl(textfile);   fprintf( ' %s\n', buffer );	% comment
    nostates = fscanf( textfile, '  nostates:\t%d\n' ); 	% number of states in polygoninference
    fprintf( '  nostates:\t\t\t%d\n', nostates );
    nomarkers = fscanf( textfile, 'markers:\t%d\n' );       % number of markers in polygoninference
    fprintf( '  nomarkers:\t\t%d\n', nomarkers );
    dim = fscanf( textfile, '  dim:\t%d\n' );               % number of dimensions in polygoninference
    fprintf( '  dim:\t\t\t\t%d\n', dim );
    nodata = fscanf( textfile, '  nodata:\t%d\n' );         % number of datapoints in polygoninference
    maxnodata = min(skip+maxnodata,nodata)-skip;            % cut short if necessary
    fprintf( '  nodata:\t\t\t%d\n', maxnodata );
    %fgetl(textfile);    % just empty lines
    %fgetl(textfile);    % just empty lines
    %fgetl(textfile);    % just empty lines
    %fgetl(textfile);    % just empty lines
    %fgetl(textfile);    % just empty lines
    % read samples:
    pd = textscan(textfile,'%f'); pd = pd{1}; nowrittendata = numel(pd)/(dim*nomarkers);
    mydata = permute( reshape( pd, [dim,nomarkers,nowrittendata] ), [1,2,3] );
    mydata = mydata(:,:,skip+(1:maxnodata));                % cut short
    nodata = size(mydata,3);
end     % end of readpolygonstates function
function[] = getfinalstatistics( truegeo, sorttype )
    % computes some basic statistics
    
    % get basic parameters:
    nochains = size(truegeo,2); if( nochains<2 ), return; end   % only run with multiple chains
    nostates = truegeo{1}.nostates; nomarkers = truegeo{1}.nomarkers;
    statsrange = truegeo{1}.statsrange; nstatsrange = numel(statsrange);
    ndatasets = size(truegeo{1}.datasetvskin,1);
    trunkfilename = truegeo{1}.subfoldername;
    withoutput = (1==1);            % for graphical output
    truegeo_sorted = sortstatelabels( truegeo, sorttype );
    
    r = -ones(nochains,nstatsrange);
    within = zeros(1,nochains);   between = zeros(1,nochains);
    for j_state = 1:nostates
        fprintf( '  ...state %d:\n',j_state );
        % lengths:
        for j_marker = 1:(nomarkers-1)
            for jj_marker = (j_marker+1):nomarkers
                minedges = +Inf; maxedges = -Inf; edgesspacing = +Inf;  % reset histogram edges for this parameter
                % ...graphical:
                fig = figure; title(sprintf('length l%d%d(%d) evolution stacked (after burnin)',j_marker,jj_marker,j_state)); xlabel('iteration');  ylabel( sprintf('l%d%d(%d)',j_marker,jj_marker,j_state) );
                hold on;
                for j_chain = 1:nochains
                    r(j_chain,:) = reshape( sqrt( sum((truegeo_sorted{j_chain}.flpos_hist(:,j_marker,j_state,statsrange)-truegeo_sorted{j_chain}.flpos_hist(:,jj_marker,j_state,statsrange)).^2,1) ), [1,nstatsrange] );
                    within(j_chain) = std(r(j_chain,:));  between(j_chain) = mean(r(j_chain,:));
                    % ....evolution:
                    plot( statsrange, r(j_chain,:) );
                    % ....histogram;
                    [~,edges_here] = histcounts( r(j_chain,:) );
                    minedges = min( minedges,edges_here(1) );   maxedges = max( maxedges,edges_here(end) );   edgesspacing = min( edgesspacing, edges_here(2)-edges_here(1) );
                end     % end of chains loop
                hold off; drawnow;  getprinted( fig, get(get(gca,'title'),'string'), trunkfilename,withoutput );
                fig = figure; title( sprintf('length l%d%d(%d) histogram stacked (after burnin)',j_marker,jj_marker,j_state) ); xlabel( sprintf('l%d%d(%d)',j_marker,jj_marker,j_state) ); ylabel('freq'); hold on;
                for j_chain = 1:nochains
                    histogram( r(j_chain,:), 'BinEdges',minedges:edgesspacing:maxedges );
                end     % end of chains loop
                hold off; drawnow;  getprinted( fig, get(get(gca,'title'),'string'), trunkfilename,withoutput );
                % ...control-window:
                GRR_simple = getGRRfrombetweenandwithin( between, within, nstatsrange );
                fprintf( '   l%d%d(%d):\t\t\t %1.5e +- %1.5e\t(GRR_simple = %1.5e)\n', j_marker,jj_marker,j_state, mean(r(:)), std(r(:)), GRR_simple );
            end 	% end of second markers loop
        end     % end of first makers loop

        % measurement errors:
        if( nomarkers>2 )
            for j_marker = 1:nomarkers
                % xy-direction:
                minedges = +Inf; maxedges = -Inf; edgesspacing = +Inf;  % reset histogram edges for this parameter
                % ...graphical:
                fig = figure; title(sprintf('meas.error merr%d_xy(%d) evolution stacked (after burnin)',j_marker,j_state), 'Interpreter','none'); xlabel('iteration');  ylabel( sprintf('merr%d_xy(%d)',j_marker,j_state), 'Interpreter','none' );
                hold on;
                for j_chain = 1:nochains
                    r(j_chain,:) = reshape( truegeo_sorted{j_chain}.merr_hist(1,1,j_marker,j_state,statsrange), [1,nstatsrange] );
                    within(j_chain) = std(r(j_chain,:));  between(j_chain) = mean(r(j_chain,:));
                    % ....evolution:
                    plot( statsrange, r(j_chain,:) );
                    % ....histogram;
                    [~,edges_here] = histcounts( r(j_chain,:) );
                    minedges = min( minedges,edges_here(1) );   maxedges = max( maxedges,edges_here(end) );   edgesspacing = min( edgesspacing, edges_here(2)-edges_here(1) );
                end     % end of chains loop
                hold off; drawnow;  getprinted( fig, get(get(gca,'title'),'string'), trunkfilename,withoutput );
                fig = figure; title( sprintf('meas.error merr%d_xy(%d) histogram stacked (after burnin)',j_marker,j_state), 'Interpreter','none' ); xlabel( sprintf('merr%d_xy(%d)',j_marker,j_state), 'Interpreter','none' ); ylabel('freq'); hold on; 
                for j_chain = 1:nochains
                    histogram( r(j_chain,:), 'BinEdges',minedges:edgesspacing:maxedges );
                end     % end of chains loop
                hold off; drawnow;  getprinted( fig, get(get(gca,'title'),'string'), trunkfilename,withoutput );
                % ...control-window:
                GRR_simple = getGRRfrombetweenandwithin( between, within, nstatsrange );
                fprintf( '   merr%d_xy(%d):\t\t %1.5e +- %1.5e\t(GRR_simple = %1.5e)\n', j_marker,j_state, mean(r(:)), std(r(:)), GRR_simple );

                % z-direction:
                minedges = +Inf; maxedges = -Inf; edgesspacing = +Inf;  % reset histogram edges for this parameter
                % ...graphical:
                fig = figure; title(sprintf('meas.error merr%d_z(%d) evolution stacked (after burnin)',j_marker,j_state), 'Interpreter','none'); xlabel('iteration');  ylabel( sprintf('merr%d_z(%d)',j_marker,j_state), 'Interpreter','none' );
                hold on;
                for j_chain = 1:nochains
                    r(j_chain,:) = reshape( truegeo_sorted{j_chain}.merr_hist(3,3,j_marker,j_state,statsrange), [1,nstatsrange] );
                    within(j_chain) = std(r(j_chain,:));  between(j_chain) = mean(r(j_chain,:));
                    % ....evolution:
                    plot( statsrange, r(j_chain,:) );
                    % ....histogram;
                    [~,edges_here] = histcounts( r(j_chain,:) );
                    minedges = min( minedges,edges_here(1) );   maxedges = max( maxedges,edges_here(end) );   edgesspacing = min( edgesspacing, edges_here(2)-edges_here(1) );
                end     % end of chains loop
                hold off; drawnow;  getprinted( fig, get(get(gca,'title'),'string'), trunkfilename,withoutput );
                fig = figure; title( sprintf('meas.error merr%d_z(%d) histogram stacked (after burnin)',j_marker,j_state), 'Interpreter','none' ); xlabel( sprintf('merr%d_z(%d)',j_marker,j_state), 'Interpreter','none' ); ylabel('freq'); hold on; 
                for j_chain = 1:nochains
                    histogram( r(j_chain,:), 'BinEdges',minedges:edgesspacing:maxedges );
                end     % end of chains loop
                hold off; drawnow;  getprinted( fig, get(get(gca,'title'),'string'), trunkfilename,withoutput );
                % ...control-window:
                GRR_simple = getGRRfrombetweenandwithin( between, within, nstatsrange );
                fprintf( '   merr%d_z (%d):\t\t %1.5e +- %1.5e\t(GRR_simple = %1.5e)\n', j_marker,j_state, mean(r(:)), std(r(:)), GRR_simple );
            end     % end of first makers loop
        else    % i.e. if only two markers
            j_marker = 1; jj_marker = 2;
            % xy-direction:
            minedges = +Inf; maxedges = -Inf; edgesspacing = +Inf;  % reset histogram edges for this parameter
            % ...graphical:
            fig = figure; title(sprintf('meas.error merr%d%d_xy(%d) evolution stacked (after burnin)',j_marker,jj_marker,j_state), 'Interpreter','none'); xlabel('iteration');  ylabel( sprintf('merr%d%d_xy(%d)',j_marker,jj_marker,j_state), 'Interpreter','none' );
            hold on;
            for j_chain = 1:nochains
                r(j_chain,:) = reshape( sqrt( truegeo_sorted{j_chain}.merr_hist(1,1,j_marker,j_state,statsrange).^2 + truegeo_sorted{j_chain}.merr_hist(1,1,jj_marker,j_state,statsrange).^2), [1,nstatsrange] );
                within(j_chain) = std(r(j_chain,:));  between(j_chain) = mean(r(j_chain,:));
                % ....evolution:
                plot( statsrange, r(j_chain,:) );
                % ....histogram;
                [~,edges_here] = histcounts( r(j_chain,:) );
                minedges = min( minedges,edges_here(1) );   maxedges = max( maxedges,edges_here(end) );   edgesspacing = min( edgesspacing, edges_here(2)-edges_here(1) );
            end     % end of chains loop
            hold off; drawnow;  getprinted( fig, get(get(gca,'title'),'string'), trunkfilename,withoutput );
            fig = figure; title( sprintf('meas.error merr%d%d_xy(%d) histogram stacked (after burnin)',j_marker,jj_marker,j_state), 'Interpreter','none' ); xlabel( sprintf('merr%d%d_xy(%d)',j_marker,jj_marker,j_state), 'Interpreter','none' ); ylabel('freq'); hold on; 
            for j_chain = 1:nochains
                histogram( r(j_chain,:), 'BinEdges',minedges:edgesspacing:maxedges );
            end     % end of chains loop
            hold off; drawnow;  getprinted( fig, get(get(gca,'title'),'string'), trunkfilename,withoutput );
            % ...control-window:
            GRR_simple = getGRRfrombetweenandwithin( between, within, nstatsrange );
            fprintf( '   merr%d%d_xy(%d):\t %1.5e +- %1.5e\t(GRR_simple = %1.5e)\n', j_marker,jj_marker,j_state, mean(r(:)), std(r(:)), GRR_simple );

            % z-direction:
            minedges = +Inf; maxedges = -Inf; edgesspacing = +Inf;  % reset histogram edges for this parameter
            % ...graphical:
            fig = figure; title(sprintf('meas.error merr%d%d_z(%d) evolution stacked (after burnin)',j_marker,jj_marker,j_state), 'Interpreter','none'); xlabel('iteration');  ylabel( sprintf('merr%d%d_z(%d)',j_marker,jj_marker,j_state), 'Interpreter','none' );
            hold on;
            for j_chain = 1:nochains
                r(j_chain,:) = reshape( sqrt( truegeo_sorted{j_chain}.merr_hist(3,3,j_marker,j_state,statsrange).^2 + truegeo_sorted{j_chain}.merr_hist(3,3,j_marker,j_state,statsrange).^2), [1,nstatsrange] );
                within(j_chain) = std(r(j_chain,:));  between(j_chain) = mean(r(j_chain,:));
                % ....evolution:
                plot( statsrange, r(j_chain,:) );
                % ....histogram;
                [~,edges_here] = histcounts( r(j_chain,:) );
                minedges = min( minedges,edges_here(1) );   maxedges = max( maxedges,edges_here(end) );   edgesspacing = min( edgesspacing, edges_here(2)-edges_here(1) );
            end     % end of chains loop
            hold off; drawnow;  getprinted( fig, get(get(gca,'title'),'string'), trunkfilename,withoutput );
            fig = figure; title( sprintf('meas.error merr%d%d_z(%d) histogram stacked (after burnin)',j_marker,jj_marker,j_state), 'Interpreter','none' ); xlabel( sprintf('merr%d%d_z(%d)',j_marker,jj_marker,j_state), 'Interpreter','none' ); ylabel('freq'); hold on; 
            for j_chain = 1:nochains
                histogram( r(j_chain,:), 'BinEdges',minedges:edgesspacing:maxedges );
            end     % end of chains loop
            hold off; drawnow;  getprinted( fig, get(get(gca,'title'),'string'), trunkfilename,withoutput );
            % ...control-window:
            GRR_simple = getGRRfrombetweenandwithin( between, within, nstatsrange );
            fprintf( '   merr%d%d_z (%d):\t %1.5e +- %1.5e\t(GRR_simple = %1.5e)\n', j_marker,jj_marker,j_state, mean(r(:)), std(r(:)), GRR_simple );
        end     % end if only two markers

        % state proportions:
        if( nostates>1 )
            for j_dset = 1:ndatasets
                minedges = +Inf; maxedges = -Inf; edgesspacing = +Inf;  % reset histogram edges for this parameter
                % ...graphical:
                fig = figure; title(sprintf('state proportion p(%d,%d) evolution stacked (after burnin)',j_state,j_dset)); xlabel('iteration');  ylabel( sprintf('p(%d,%d)',j_state,j_dset) );
                hold on;
                for j_chain = 1:nochains
                    r(j_chain,:) = truegeo_sorted{j_chain}.statefrac_hist(j_state,j_dset,statsrange);
                    within(j_chain) = std(r(j_chain,:));  between(j_chain) = mean(r(j_chain,:));
                    % ....evolution:
                    plot( statsrange, r(j_chain,:) );
                    % ....histogram;
                    [~,edges_here] = histcounts( r(j_chain,:) );
                    minedges = min( minedges,edges_here(1) );   maxedges = max( maxedges,edges_here(end) );   edgesspacing = min( edgesspacing, edges_here(2)-edges_here(1) );
                end     % end of chains loop
                hold off; drawnow;  getprinted( fig, get(get(gca,'title'),'string'), trunkfilename,withoutput );
                fig = figure; title( sprintf('state proportion p(%d,%d) histogram stacked (after burnin)',j_state,j_dset) ); xlabel( sprintf('p(%d,%d)',j_state,j_dset) ); ylabel('freq'); hold on; 
                for j_chain = 1:nochains
                    histogram( r(j_chain,:), 'BinEdges',minedges:edgesspacing:maxedges );
                end     % end of chains loop
                hold off; drawnow;  getprinted( fig, get(get(gca,'title'),'string'), trunkfilename,withoutput );
                % ...control-window:
                GRR_simple = getGRRfrombetweenandwithin( between, within, nstatsrange );
                fprintf( '   p(%d,%d):\t\t\t %1.5e +- %1.5e\t(GRR_simple = %1.5e)\n', j_state,j_dset, mean(r(:)), std(r(:)), GRR_simple );
            end     % end of datasets loop
        end     % end if multiple states
    end     % end of states loop
end     % end of getfinalstatistics function
function[ truegeo_sorted ] = sortstatelabels( truegeo, sorttype )
    % sorts the statelabels according to some criteria
    
    % get basic parameters:
    nochains = size(truegeo,2); if( nochains<2 ), truegeo_sorted = truegeo; return; end               % only run with multiple chains
    fprintf( ' Info - sortstatelabels: Type is %s.\n', sorttype );
    if( strcmp(sorttype,'none') )
        truegeo_sorted = truegeo; return;
    elseif( strcmp(sorttype,'circumference') )
        mysorttype = 1;                         % internal codification
    elseif( strcmp(sorttype,'proportion') )
        mysorttype = 2;                         % internal codification
    else    % unknown 
        fprintf( ' Info - sortstatelabels: Type is not known - abort sorting.\n' ); truegeo_sorted = truegeo; return;
    end     % end of going through various sorttypes
    nostates = truegeo{1}.nostates; nomarkers = truegeo{1}.nomarkers;       % should be same for all independent chains
    MCmax = size(truegeo{1}.flpos_hist,4);      % should be same for all independent chains
    truegeo_sorted{nochains} = [];              % allocate memory
    comparable{nostates} = [];                  % groups comparable states together
    
    % get comparable:
    j_comp = 0;                                 % initialise counter
    while( j_comp<size(comparable,2) )          % stop, if all states are affiliated
        j_comp = j_comp + 1;                    % go through next entry
        comparable{j_comp} = j_comp;            % will get removed afterwards, if redundant
        for jj_comp = 1:(j_comp-1)
            cond = all( truegeo{1}.posstatevskin(j_comp,:)==truegeo{1}.posstatevskin(jj_comp,:), 'all' );
            cond = cond & all( truegeo{1}.lengthpriordensity.pars{j_comp}==truegeo{1}.lengthpriordensity.pars{jj_comp}, 'all' );
            cond = cond & all( truegeo{1}.shapepar_prior{j_comp}==truegeo{1}.shapepar_prior{jj_comp}, 'all' );
            cond = cond & all( truegeo{1}.merrmeans_prior{j_comp}==truegeo{1}.merrmeans_prior{jj_comp}, 'all' );
            if( cond )
                comparable{jj_comp} = [comparable{jj_comp},j_comp]; % add to list of states with same prior
                comparable(j_comp) = []; j_comp = j_comp-1; % remove as an independent entry
                break;                        	% don't have to go through rest of state anymore
            end     % end if same prior condition
        end     % end of going through subloop
    end     % end of determining comparable
    fprintf( ' Info - sortstatelabels: Comparable priors for states:\t' );
    for j_comp = 1:size(comparable,2)
        fprintf( '[ %s]', sprintf('%d ',comparable{j_comp}) );
    end     % end of comparable loop
    fprintf('\n');
    
    % sort within same prior:
    for j_comp = 1:size(comparable,2)
        for j_chain = 1:nochains
            for j_MC = 1:MCmax
                % ...get crit for this iteration:
                crit = zeros(size(comparable{j_comp}));     % reset; row-vector
                for j_crit = 1:numel(crit)
                    j_state = comparable{j_comp}(j_crit);   % short-hand notation
                    if( mysorttype==1 )
                        for j_marker = 1:(nomarkers-1)
                            for jj_marker = (j_marker+1):nomarkers
                                crit(j_crit) = crit(j_crit) + squeeze(sqrt( sum((truegeo{j_chain}.flpos_hist(:,j_marker,j_state,j_MC)-truegeo{j_chain}.flpos_hist(:,jj_marker,j_state,j_MC)).^2,1) ));
                            end     % end of inner markers loop
                        end     % end of outer markers loop
                    elseif( mysorttype==2 )
                        crit(j_crit) = truegeo{j_chain}.statefrac_hist(j_state,1,j_MC);
                    end     % end of distinguishing mysorttype
                end     % end of getting crit
                % ...get ordering for this iteration:
                [~,ordering] = sort( crit ); newstateorder = comparable{j_comp}(ordering);
                % ...get trugeo_sorted for this iteration:
                truegeo_sorted{j_chain}.flpos_hist(:,:,comparable{j_comp},j_MC) = truegeo{j_chain}.flpos_hist(:,:,newstateorder,j_MC);
                truegeo_sorted{j_chain}.merr_hist(:,:,:,comparable{j_comp},j_MC) = truegeo{j_chain}.merr_hist(:,:,:,newstateorder,j_MC);
                truegeo_sorted{j_chain}.statefrac_hist(comparable{j_comp},:,j_MC) = truegeo{j_chain}.statefrac_hist(newstateorder,:,j_MC);
                truegeo_sorted{j_chain}.statevskin_hist(comparable{j_comp},:,j_MC) = truegeo{j_chain}.statevskin_hist(newstateorder,:,j_MC);
            end     % end of MCit loop
        end     % end of chains loop
    end     % end of comparable loop
end     % end of sortstatelabels function
function[] = getprinted( fig, pictitle, subfoldername,withoutput )
    % exports figures in various formats

    if( withoutput )
        if( iscell(pictitle)==1 ), pictitle = pictitle{1}; end
        figfilename = fullfile( subfoldername, [pictitle,'.png'] );
        print( fig, figfilename, '-dpng' );
        figfilename = fullfile( subfoldername, [pictitle,'.fig'] );
        savefig( fig, figfilename );

        % close picture again after a short while:
        pause(0.3);
        close(fig);
    end     % end if withfileexport
end     % end of getprinted function
function[ GRR_simple ] = getGRRfrombetweenandwithin( between, within, nstatsrange )
    
    W = mean(within.^2);    B = nstatsrange*var(between);
    V_simple = (1-1/nstatsrange)*W + (1/nstatsrange)*B;
    
    % get output:
    GRR_simple = sqrt( V_simple./W );
end     % end of getGRRfrombetweenandwithin function
