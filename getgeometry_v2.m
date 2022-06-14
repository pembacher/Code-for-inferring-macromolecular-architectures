classdef getgeometry_v2 < handle
    % gives underlying geometry for a set of samples of GPS sister pairs
    % assumes equal geomtry of each sister and only twist around fl1
    % (in all three dimensions);
    % also assumes different measurements only differ by different
    % perspectives (individual translations, rotations) but otherwise
    % identical geometries
    
    properties ( Access=public )
        % physical/experimental data:
        inputfilename = [];     % inputfilename
        flpos_data = [];        % observed positions
        posstatevskin = -ones(1,1); % nostates x nodata sized matrix of logicals, assigning which state to be present in which measurement
        datasetvskin = -ones(1,1);  % nodatasets x nodata sized matrix of logicals, assigning which measurement to be in which dataset
        % initial settings:
        flpos_init = zeros(3,1);% starting position of MCMC for fluorophore positions
        merr_init = zeros(3,3); % starting position of MCMC for measurement errors
        lengthpriortype = 1;    % '0' for flat prior on positions, '1' for flat prior on lengths
        lengthpriordensity = [];% defines density with respect to flat prior on lengths, if '1' is selected in lengthpriortype
        merrmeans_prior = [];   % prior in nm for means of Gamma-distribution of merr; 0 is flat
        shapepar_prior = [];    % prior on shapeparameter of gamma distribution; 1-3/2 is flat
        
        % numerical data:
        globaltimestamp = [];	% clock for global reference
        burnin = 1;             % first couple of MC-iterations that gets ignored
        MCmax = 1;              % maximal number of MCMC steps
        subsamplefreq = 1;      % frequency with which samples are added to the history
        stepsizes = ones(11,1); % list of stepsizes in order they are used for updateS
        errorupdatetype = 'Gibbs';  % 'rwMH' or 'Gibbs'
        nostates = 1;           % indicates number of states that each measurement can be in
        nodsets = 1;            % indicates number of datasets that each measurement can be in 
        outputfreq = 1000;      % frequency with which current samples are output in the controlwindow
        comment = '';           % comment in output
    end     % end of entirely public properties
    properties( SetAccess = protected, GetAccess = public ) % i.e. protected to write; public to read
        % physcial parameters:
        flpos_curr = [];                % current sample of template position
        merr_curr = [];                 % current sample of measurement error
        trans_curr = zeros(3,1);        % current translations
        rot_curr = zeros(3,3,1);        % current rotations
        flpos_pred_curr = [];           % current true position
        flpos_prop = [];                % proposal of template position
        merr_prop = [];                 % proposal of measurement error
        trans_prop = zeros(3,1);        % proposal of translations
        rot_prop = zeros(3,3,1);        % proposal of rotations
        flpos_pred_prop = {};           % proposal of true position
        logtarget_curr = 0;     logtarget_prop = 0;     % logarithm of current/proposed target value; gives sarw-rescaled density unlike auxpars_sarw.logtarget_curr; needed for mixture model of multistate
        statefrac_curr = 0;     statefrac_prop = 0;     % proportion of how often a state shows up
        statevskin_curr = 0;    statevskin_prop = 0;    % gives which state each kinetochore is in (right sister); notrianglestatesxnodata sized matrix
        % auxiliary parameters:
        nomarkers = 1;              % number of fluorophores
        nodata = 1;                 % number of experiments (i.e. length of each fl-vector)
        MCit = 0;                   % current iteration
        t1 = 0;                     % time, when values are initialised
        nextoptstepstime = 1;       % next time stepsizes are optimised
       	rejected = zeros(11,1);   	% number of rejected samples for [true geometry, measurement error, perspectives,large-scale-rotations_Gibbs, persp_prior_phik, persp_prior_phiaxis, alljoint, indep, statechange, chromshift]
        samplecounter = zeros(11,1);% number of total samples for [true geometry, measurement error, perspectives,large-scale-rotations_Gibbs, persp_prior_phik, persp_prior_phiaxis, alljoint, indep, statechange, chromshift]   (only gets updated for rw oae samplers)
        adjustedtoofar = zeros(11,1);% memory how often optimal rejection-rate have been passed-over
        subfoldername = 'mynewfolder';  % folder for output
        statsrange = [];            % range of indices in hist-variables, that are evaluated as post-burnin
        without = 1;                % '0' for no output, except warnings; '1' for basic output; '2' for mode details etc
        
        % statistical parameters:
        flpos_hist = [];            % for template positions
        merr_hist = [];             % for measurement errors
        trans_hist = zeros(3,1);    % for translations
        rot_hist = zeros(3,3,1);    % for rotations
        statefrac_hist = zeros(1,1,1);  % for state proportions
        statevskin_hist = zeros(1,1,1); % first index is state-index, second index is data-index, third index is MCit index
        logtarget_hist = [];        % for logtarget values
    end     % end of public-read properties
    
    methods( Access=public )
        function[] = initialiseparameters( obj )
            obj.nomarkers = size( obj.flpos_data,2 );
            obj.nodata = size( obj.flpos_data,3 );
            
            % initialise length prior:  (this is probably wrong prior for external prior, so need to disable priors when updating draw from external prior)
            if( (obj.lengthpriortype==1) && (isempty(obj.lengthpriordensity)) ) % i.e. prior on lengths and no density given manually
                obj.lengthpriordensity.type = 1;    % '1' for intervals
                for j_state = 1:obj.nostates
                    obj.lengthpriordensity.pars{j_state} = [ zeros(obj.nomarkers*(obj.nomarkers-1)/2,1), Inf(obj.nomarkers*(obj.nomarkers-1)/2,1) ];    % columns go through markers in order AB, AC, ..., BC, ...
                end     % end of states loop
            end     % end of giving default prior on lengths
            
            % initialise error priors:
            if( isempty(obj.shapepar_prior) )	% i.e. no prior given at all
                for j_state = 1:obj.nostates
                    obj.shapepar_prior{j_state} = ones(obj.nomarkers,1)*(1-(3/2));	% shapeparameter of prior of Gamma-distribution of merr;       1-3/2 is flat; column vector for each marker
                end     % end of states loop
            end     % end of setting shapepar prior
            if( isempty(obj.merrmeans_prior) )% i.e. no prior given at all
                for j_state = 1:obj.nostates
                    obj.merrmeans_prior{j_state} = ones(obj.nomarkers,1)*0*[1,1,2];	% prior in nm for means of Gamma-distribution of merr; 0 is flat; column-index for each marker, row-index for the three space dimensions
                end     % end of states loop
            end     % end of setting merrmeans prior
            
            % initialises the position- and error-parameters:
            bestguess = (0==1);
            obj.flpos_curr = zeros(3,obj.nomarkers,obj.nostates);
            obj.merr_curr = zeros(3,3,obj.nomarkers,obj.nostates);
            validarchitecture = (0==1);             % becomes false, if no physically possible geometry is sampled
            while( ~validarchitecture )
                for j_state = 1:obj.nostates
                    if( bestguess==1 )                  % use best guess to initialise
                        obj.flpos_curr(:,:,j_state) = obj.flpos_init(:,:,j_state);% current sample of template position
                        obj.merr_curr(:,:,j_state) = obj.merr_init(:,:,j_state);  % current sample of measurement error
                    else                                % use prior to initialise
                        if( obj.lengthpriortype==0 )                % i.e. prior on positions
                            obj.flpos_curr(:,1,j_state) = zeros(3,1);% current sample of template position of fl1
                            obj.flpos_curr(:,2:end,j_state) = 300*rand([3,obj.nomarkers-1]);% current sample of true position of fl2
                        elseif( obj.lengthpriortype==1 )    % prior on lengths; choose uniformly in prior
                            j_length = 0;                       % goes through length indices in order
                            for j_marker = 1:obj.nomarkers
                                if( j_marker==1 )
                                    obj.flpos_curr(:,j_marker,j_state) = zeros(3,1);    % in origin
                                elseif( j_marker==2 )
                                    j_length = j_length + 1;	% need to evaluate one more length
                                    if( all(~isinf(obj.lengthpriordensity.pars{j_state}(j_length,:))) )	% i.e. no infs
                                        r12 = obj.lengthpriordensity.pars{j_state}(j_length,1) + (2*rand()-1)*obj.lengthpriordensity.pars{j_state}(j_length,2);
                                    else
                                        r12 = 300*rand();
                                    end     % end if length prior
                                    obj.flpos_curr(:,j_marker,j_state) = [1;0;0]*r12;   % on positive x-axis
                                elseif( j_marker==3 )
                                    j_length = j_length + 1;	% need to evaluate one more length
                                    if( all(~isinf(obj.lengthpriordensity.pars{j_state}(j_length,:))) )	% i.e. no infs
                                        r13 = obj.lengthpriordensity.pars{j_state}(j_length,1) + (2*rand()-1)*obj.lengthpriordensity.pars{j_state}(j_length,2);
                                    else
                                        r13 = 300*rand();
                                    end     % end if length prior
                                    j_length = j_length + 1;	% need to evaluate one more length
                                    if( all(~isinf(obj.lengthpriordensity.pars{j_state}(j_length,:))) )	% i.e. no infs
                                        r23 = obj.lengthpriordensity.pars{j_state}(j_length,1) + (2*rand()-1)*obj.lengthpriordensity.pars{j_state}(j_length,2);
                                    else
                                        r23 = 300*rand();
                                    end     % end if length prior
                                    cosalpha = (r12^2+r13^2-r23^2)/(2*r12*r13);     % angle 213
                                    obj.flpos_curr(:,j_marker,j_state) = [ cosalpha*r13; sin(abs(acos(cosalpha)))*r13; 0 ];	% abs, so always in upper half plane of xy-plane
                                elseif( j_marker>3 )
                                    obj.flpos_curr(:,j_marker,j_state) = 300*rand(3,1); % random for now
                                end     % end of going through possible markers
                            end     % end of markers loop
                        end     % end if lengthpriortype
                        for j_marker = 1:obj.nomarkers
                            if( obj.shapepar_prior{j_state}(j_marker)<0 )
                                obj.merr_curr(:,:,j_marker,j_state) = diag(rand(3,1))*200;	% current sample of measurement error of fl1(s)
                            else
                                obj.merr_curr(:,:,j_marker,j_state) = diag( 1./sqrt( gamrnd( obj.shapepar_prior{j_state}(j_marker), 1./(obj.merrmeans_prior{j_state}(j_marker,:)*sqrt(obj.shapepar_prior{j_state}(j_marker))).^2 ) ) );
                            end
                            obj.merr_curr(2,2,j_marker,j_state) = obj.merr_curr(1,1,j_marker,j_state);    % ignore second component and identify with first one
                        end     % end of markers loop
                    end     % end if bestguess
                end     % end of state-loop
                if( (~all(~isinf(obj.flpos_curr),'all')) || (~all(~isnan(obj.flpos_curr),'all')) || (~all(isreal(obj.flpos_curr),'all')) )
                    validarchitecture = (0==1);     % something went wrong
                else
                    validarchitecture = (~isinf(obj.getabsoluteloglengthprior( obj.flpos_curr ))); % try again, if pathological for this prior
                end     % end if pathological
            end     % end if validarchitecture
            % copy over to inits:
            obj.flpos_init = obj.flpos_curr;
            obj.merr_init = obj.merr_curr;
                        
            % initialise perspectives:
            obj.trans_curr = zeros(3,obj.nodata);
            obj.rot_curr = zeros(3,3,obj.nodata);
            for j_data = 1:obj.nodata
                obj.trans_curr(:,j_data) = 2*(rand(3,1)-0.5)*2e4;   % in large box
                obj.rot_curr(:,:,j_data) = obj.sampleisorotation();
            end     % end going through data
            
            % initialise triangle state labels:
            obj.statefrac_curr = rand(obj.nostates,obj.nodsets); obj.statefrac_curr = obj.statefrac_curr./sum(obj.statefrac_curr,1);
            
            % initialise triangle state labels:
            obj.statevskin_curr = false(obj.nostates,obj.nodata);
            while( min(sum(obj.statevskin_curr,2))<3 )  % hard boundary for state occupations
                cindex = obj.nostates*(0:(obj.nodata-1));
                for j_dset = 1:obj.nodsets
                    data_here = obj.datasetvskin(j_dset,:); nodata_here = sum(data_here);
                    cindex_here = cindex(data_here);
                    cumsumfracs = repmat([0;obj.statefrac_curr(:,j_dset)],[1,nodata_here]);
                    cumsumfracs = cumsumfracs.*[zeros(1,nodata_here);obj.posstatevskin(1:obj.nostates,data_here)]; 
                    cumsumfracs = cumsum(cumsumfracs,1); 
                    cumsumfracs = cumsumfracs./cumsumfracs(end,:);   % cumsum of all possible states (impossible ones are set to zero probability)
                    stateoptions = repmat((1:(obj.nostates+1))',[1,nodata_here]);
                    buffer = (rand(1,nodata_here)>=cumsumfracs); buffer = max(stateoptions.*buffer,[],1); obj.statevskin_curr(cindex_here(buffer>0)+buffer(buffer>0)) = (1==1);    % note: buffer==0, if no state is possible (i.e. posstatevskin==0 for all states), so want to keep (0==1) entry there
                end     % end of datasets loop
            end     % end keep trying until at least three measurements in each state
            
            % take current samples as proposals:
            obj.flpos_prop = obj.flpos_curr;
            obj.merr_prop = obj.merr_curr;
            obj.trans_prop = obj.trans_curr;
         	obj.rot_prop = obj.rot_curr;
         	obj.statefrac_prop = obj.statefrac_curr;
            obj.statevskin_prop = obj.statevskin_curr;
            
            % compute current predictions:
            % ...first compute proposed ones (due to notation of updates):
            obj.flpos_pred_prop = zeros(3,obj.nomarkers,obj.nodata);
            obj.updatepositionpredictions(-1,-1);
            % ...and copy into current values:
            obj.flpos_pred_curr = obj.flpos_pred_prop;
            
            
            % set auxiliary parameters:
            obj.MCit = 0;                   % initial value counts as zero
            noups = obj.nostates*obj.nomarkers + strcmp(obj.errorupdatetype,'rwMH')*obj.nostates*obj.nomarkers + 1 + (obj.nostates>1);
            obj.rejected = zeros(noups,1);	% no rejections yet
            obj.samplecounter = zeros(noups,1); % no samples yet
            obj.stepsizes = ones(noups,1);  % stepsizes
            obj.adjustedtoofar = zeros(noups,1);% memory how often optimal rejection-rate have been passed-over
            obj.t1 = now;                   % current time at end of initialisation
            if( isempty(obj.globaltimestamp) ), obj.globaltimestamp=clock(); end
            if( 1==1 )
                obj.subfoldername = sprintf( '%d-%02.0f-%02.0f_%02.0f-%02.0f_%s', obj.globaltimestamp(1:5), 'getgeometry output' );
                if( ~exist(obj.subfoldername,'dir') )
                    mkdir(obj.subfoldername);               % creates a subfolder of the given name in the current directory
                end     % end if directory does not already exist
            end     % end if with fileexport
            
            % compute current target value:
            obj.logtarget_prop = obj.getabsolutetargetprobability(-1,-1);
            obj.logtarget_curr = obj.logtarget_prop;
            
            % set statistic parameters:
            obj.statsrange = (obj.burnin + 1):obj.MCmax;
            for j_state = 1:obj.nostates
                obj.flpos_hist = zeros(3,obj.nomarkers,obj.nostates,obj.MCmax);    	obj.flpos_hist(:,:,:,1) = obj.flpos_curr;
                obj.merr_hist = zeros(3,3,obj.nomarkers,obj.nostates,obj.MCmax);    obj.merr_hist(:,:,:,:,1) = obj.merr_curr;
            end     % end of going through states
            obj.trans_hist = zeros(3,obj.nodata,obj.MCmax);                         obj.trans_hist(:,:,1) = obj.trans_curr;
            obj.rot_hist = zeros(3,3,obj.nodata,obj.MCmax);                         obj.rot_hist(:,:,:,1) = obj.rot_curr;
            obj.statefrac_hist = zeros(obj.nostates,obj.nodsets,obj.MCmax);        	obj.statefrac_hist(:,:,1) = obj.statefrac_curr;
            obj.statevskin_hist = false(obj.nostates,obj.nodata,obj.MCmax);         obj.statevskin_hist(:,:,1) = obj.statevskin_curr;
            obj.logtarget_hist = zeros(1,obj.MCmax);
        end     % end of initialiseparameters function
        function[] = control( obj, runrange )
            % runs the MCMC
            
            %diary 'getgeometry controlwindow output.txt'; diary on; % output the following to the controlwindow
            if( nargin<2 )                              % 'obj' does count towards nargin
                runrange = 1:obj.MCmax;
                fprintf( ' (%s) Info - control (%d): Set runrange myself ([%d..%d]).\n', obj.comment,obj.MCit, runrange(1),runrange(end) );
            end     % end if runrange not given
            
            % initial output:
            if( runrange(1)==1 )                        % runrange just starts
                % ...in control-window:
                fprintf( ' (%s) Start getgeometry_v2:\n', obj.comment );
                fprintf( ' (%s) timestamp:\t\t\t\t%4.4d_%2.2d_%2.2d-%2.2d_%2.2d\n', obj.comment, obj.globaltimestamp(1:5) );
                fprintf( ' (%s) statsrange:\t\t\t[%d..%d]\n', obj.comment, obj.statsrange(1),obj.statsrange(end) );
                fprintf( ' (%s) subsamplefreq:\t\t\t%d\n', obj.comment, obj.subsamplefreq );
                fprintf( ' (%s) number of states:\t\t%d\n', obj.comment, obj.nostates );
                fprintf( ' (%s) lengthpriortype:\t\t%d\n', obj.comment, obj.lengthpriortype )
                if( obj.lengthpriortype==1 )	% i.e. prior on lengths
                    fprintf(  ' (%s) ...densitytype:\t\t%d\n', obj.comment, obj.lengthpriordensity.type );
                    for j_state = 1:obj.nostates
                        fprintf(  ' (%s) ...state %d:\t\t\t[ %s]\n', obj.comment, j_state, sprintf( '%1.3e+-%1.3e ',reshape(obj.lengthpriordensity.pars{j_state}',[1,2*obj.nomarkers*(obj.nomarkers-1)/2])) );
                    end     % end of states loop
                end     % end of outputing density parameters
                fprintf( ' (%s) errorupdatetype:\t\t%s\n', obj.comment, obj.errorupdatetype );
                fprintf( ' (%s) errorupdatepior:\n', obj.comment );
                for j_state = 1:obj.nostates
                    fprintf( ' (%s) ...state %d:\t\t\t', obj.comment, j_state );
                    for j_marker = 1:obj.nomarkers
                        fprintf( 'Gamma(%+1.3e,[%1.3e,%1.3e,%1.3e]) ', obj.shapepar_prior{j_state}(j_marker),obj.merrmeans_prior{j_state}(j_marker,:) );
                    end     % end of markers loop
                    fprintf('\n');
                end     % end of states loop
                fprintf( ' (%s) possible states for each measurement:\n', obj.comment );
                for j_state = 1:size(obj.posstatevskin,1)
                    fprintf( '     %s\n', sprintf('%d ',obj.posstatevskin(j_state,:)) );
                end     % end of states loop
                fprintf( ' (%s) dataset for each measurement:\n', obj.comment );
                for j_dset = 1:size(obj.datasetvskin,1)
                    fprintf( '     %s\n', sprintf('%d ',obj.datasetvskin(j_dset,:)) );
                end     % end of states loop
                fprintf( ' (%s) Info - control (%d): Initial state:\n', obj.comment,obj.MCit );
                obj.regularcontrolwindowoutput();
            end     % end checking if at beginning
            
            % run MCMC:
            for currentit = runrange
                obj.MCit = currentit;
                if( obj.MCit==obj.burnin+1 )
                    fprintf( ' (%s) optimised stepsizes (%d): [ %s]\n', obj.comment,obj.MCit+1, sprintf('%+6.4e ',obj.stepsizes) );
                end     % end of checking end of burnin
                obj.updatecurrentsample();
                obj.optimisestepsizes();        % must be placed after updatecurrentsample
                obj.regularcontrolwindowoutput();
            end     % end of MCMC updates
            
            % final output:
            if( runrange(end)==obj.MCmax )
                fprintf( ' (%s) Final statistics (%d):\n', obj.comment,obj.MCit );
                if( 1==1 )  % save structure externally
                    save( sprintf('%s/%d-%02.0f-%02.0f_%02.0f-%02.0f_getgeom_%s', obj.subfoldername,obj.globaltimestamp(1:5),obj.comment), 'obj' );
                end
                fprintf( ' (%s) stepsizes: [ %s]\n', obj.comment, sprintf('%+6.4e ',obj.stepsizes) );
                fprintf( ' (%s) rej rates: [ %s]\n', obj.comment, sprintf('%+6.4e ',obj.rejected./obj.samplecounter) );
                obj.outputlengthstatistics( 1, 0, 0, 0 );           % arguments are 'withcw, withhist, withscat, withevol' in this order
                obj.outputinternalanglesstatistics( 1, 0, 0, 0 );	% arguments are 'withcw, withhist, withscat, withevol' in this order
                obj.outputmarginalstates( 1, 0 );                   % arguments are 'withcw, withhist' in this order
                fprintf( ' (%s) Total time consumption of getgeometry_v2:\t\t%1.3f sec\n', obj.comment, (now-obj.t1)*24*3600 );
                %diary 'getgeometry controlwindow output.txt'; diary off;
            end     % end of if end of runrange
        end     % end of control function
        function[] = getprinted( obj, fig, pictitle )
            % exports figures in various formats
            
            if( 1==1 )
                if( iscell(pictitle)==1 ), pictitle = pictitle{1}; end
                figfilename = fullfile( obj.subfoldername, [pictitle,'_',obj.comment,'.png'] );
                print( fig, figfilename, '-dpng' );
                figfilename = fullfile( obj.subfoldername, [pictitle,'_',obj.comment,'.fig'] );
                savefig( fig, figfilename );
                
                % close picture again after a short while:
                pause(0.3);
                close(fig);
            end     % end if withfileexport
        end     % end of getprinted function
        function[] = copypropstocurr( obj )
            % copies the currently proposed values to the current values
            
            obj.flpos_curr = obj.flpos_prop;
            obj.merr_curr = obj.merr_prop;
            obj.trans_curr = obj.trans_prop;
            obj.rot_curr = obj.rot_prop;
            obj.flpos_pred_curr = obj.flpos_pred_prop;
            obj.logtarget_curr = obj.logtarget_prop;
            obj.statefrac_curr = obj.statefrac_prop;
            obj.statevskin_curr = obj.statevskin_prop;
        end     % end of copypropstocurr function
        function[] = copycurrtoprops( obj )
            % copies the current values to the currently proposed values
            
            obj.flpos_prop = obj.flpos_curr;
            obj.merr_prop = obj.merr_curr;
            obj.trans_prop = obj.trans_curr;
            obj.rot_prop = obj.rot_curr;
            obj.flpos_pred_prop = obj.flpos_pred_curr;
            obj.logtarget_prop = obj.logtarget_curr;
            obj.statefrac_prop = obj.statefrac_curr;
            obj.statevskin_prop = obj.statevskin_curr;
        end     % end of copycurrtoprops function
        function[] = changeprops( obj, flpos_in, merr_in, trans_in,rot_in, statefrac_in,statevskin_in, logtarget_in )
            obj.flpos_prop = flpos_in;
            obj.merr_prop = merr_in;
            obj.trans_prop = trans_in;
            obj.rot_prop = rot_in;
            obj.statefrac_prop = statefrac_in;
            obj.statevskin_prop = statevskin_in;
            obj.logtarget_prop = logtarget_in;
        end     % end of changeprops function
        function[ logtargetratio,logtargetprop ] = getgettargetprobabilityratio( obj, j_data, j_state )
            [logtargetratio,logtargetprop] = obj.gettargetprobabilityratio( j_data,j_state );
        end     % end of getgettargetprobabilityratio function
        function[ target ] = getgetabsolutetarget( obj, j_data, j_state )
            target = obj.getabsolutetargetprobability( j_data,j_state );
        end     % end of getgetabsolutetarget function
        function[] = getupdatepositionpredictions( obj, j_data, j_state )
            obj.updatepositionpredictions( j_data, j_state );
        end     % end of getupdatepositionpredictions function
        function[ flpos_pred_hist ] = getpredictionhist( obj, dataindex )
            % gives a history for all predictions
            
            % set auxiliary patameters:
            if( obj.MCit<obj.burnin+1 ) % still during burnin
                statsrangehere = 1:obj.MCit;         % don't include latest one, just in case it's before the update
            else                        % only considre post-burnin data
                statsrangehere = obj.statsrange(1):obj.MCit;
            end
            if( dataindex<0 )           % must be single index
                datarange = 1:obj.nodata;
            else
                datarange = dataindex;
            end
            flpos_pred_hist = zeros(3,numel(datarange),numel(statsrangehere));
            stateoptions = 1:obj.nostates;  % possible trianglestates
            
            for j_data = 1:numel(datarange)
                jj_data = datarange(j_data);
                for j_MCit = 1:numel(statsrangehere)
                    jj_MCit = statsrangehere(j_MCit);
                    if( obj.nostates>1 )
                        j_state = stateoptions(obj.statevskin_hist(:,jj_data,jj_MCit));
                    else
                        j_state = 1;
                    end     % end if multiple triangle states in mixture
                    if( ~isempty(j_state) )   % i.e. there exists a state for this datapoint
                        obj.flpos_prop(:,:,j_state) = obj.flpos_hist(:,:,j_state,jj_MCit);
                        obj.trans_prop(:,jj_data) = obj.trans_hist(:,jj_data,jj_MCit);
                        obj.rot_prop(:,:,jj_data) = obj.rot_hist(:,:,jj_data,jj_MCit);
                        obj.updatepositionpredictions(jj_data,-1);
                        flpos_pred_hist(:,:,j_data,j_MCit) = obj.flpos_pred_prop(:,:,jj_data);
                    else                                            % i.e. there is not state for this datapoint
                        flpos_pred_hist(:,j_data,j_MCit) = NaN(3,1);
                    end     % end if there exists a state for this datapoint
                end     % end of going through iterations
            end     % end of going through datapoints
            
        end     % end of getpredictionhist function
    end     % end of entirely public methods
    methods( Access=private )
        function[] = updatecurrentsample( obj )
            % gets proposals and accepts/rejects them
            
            if( obj.without>=2 )
                fprintf( ' %d, %s: %1.10e, %1.10e\n', obj.MCit, obj.comment, obj.logtarget_curr, obj.getabsolutetargetprobability(-1,-1) );
            end
            if( obj.without>=4 )
            	obj.outputcurrentandproposedpars();
            end
            %fprintf( '___update: MCit = %d____\n', obj.MCit );
            
            for j_subsample = 1:obj.subsamplefreq
                obj.copycurrtoprops();      % just to make sure sequential updates are actually sequential without remainders from the past
                j_up = 0;   % update counter (for stepsizes, rejected, samplecounter)
                % template position update:
                if( obj.without>=2 )
                    fprintf( ' Info - updatecurrentsample (%s)(%d): Positions\n', obj.comment,obj.MCit );
                    obj.comparecurrentandproposedpars();
                end
                for j_state = 1:obj.nostates
                    for j_marker = 1:obj.nomarkers                     	% go through each set of to-be-updated truegeomtry parameter(s) (see updatetruegeometryprops)
                        j_up = j_up + 1;
                        obj.updatetruegeometryprops(j_marker,j_state, j_up);
                        obj.updatepositionpredictions(-1,j_state);    	% update predicted positions, due to new true geometry
                        [logtargetratio,obj.logtarget_prop] = obj.gettargetprobabilityratio(-1,j_state);   % compute target probability of proposal  
                        if( rand()<exp(logtargetratio) )              	% accept
                            if( obj.without>=3 )
                                fprintf( ' Info - updatecurrentsample (%d)(%s): accept position (st=%d)(fl=%d)(%1.3e)\n', obj.MCit, obj.comment, j_state,j_marker, logtargetratio );
                            end
                            % ...copy proposals into current samples:
                            obj.copypropstocurr();
                        else                                  	% reject
                            if( obj.without>=3 )
                                fprintf( ' Info - updatecurrentsample (%d)(%s): reject position (st=%d)(fl=%d)(%1.3e)\n', obj.MCit, obj.comment, j_state,j_marker, logtargetratio );
                            end
                            obj.copycurrtoprops();              % needed to avoid interference with future markers
                            % no need to return to previous proposals/current values, as this is done at beginning of next sequential step or error-update, respectively
                            obj.rejected(j_up) = obj.rejected(j_up) + 1/(obj.nomarkers*obj.nostates);	% one more rejected proposal
                        end     % end of checking acceptance
                        obj.samplecounter(j_up) = obj.samplecounter(j_up) + 1/(obj.nomarkers*obj.nostates);
                    end     % end of going through all truegeometry parameters
                end     % end of state loop
                obj.standardisetruegeometry(); obj.copycurrtoprops();   % standardisetruegeometry only acts on current parameters
                
                % check new measurement error proposal:
                if( obj.without>=2 )
                    fprintf( ' Info - updatecurrentsample (%s)(%d): Error\n', obj.comment,obj.MCit );
                    obj.comparecurrentandproposedpars();
                end
                if( strcmp(obj.errorupdatetype,'rwMH') )
                    for j_state = 1:obj.nostates    % go through each state separately
                        for j_marker = 1:obj.nomarkers  % go through all fluorophore markers
                            j_up = j_up + 1;
                            obj.updateerrorprop( j_marker, j_state, j_up );
                            [logtargetratio,obj.logtarget_prop] = obj.gettargetprobabilityratio(-1,j_state);	% compute target probability of proposal
                            if( rand()<exp(logtargetratio) )                                            % accept
                                if( obj.without>=3 )
                                    fprintf( ' Info - updatecurrentsample (%d)(%s): accept error (st=%d)(fl=%d) (logtargetratio = %1.5e)\n', obj.MCit,obj.comment, j_state,j_marker, logtargetratio );
                                end
                                obj.merr_curr(:,:,j_marker,j_state) = obj.merr_prop(:,:,j_marker,j_state);
                                obj.copypropstocurr();
                            else                                                                % reject
                                if( obj.without>=3 )
                                    fprintf( ' Info - updatecurrentsample (%d)(%s): reject error (st=%d)(fl=%d) (logtargetratio = %1.5e)\n', obj.MCit,obj.comment, j_state,j_marker, logtargetratio );
                                end
                                obj.merr_prop(:,:,j_marker,j_state) = obj.merr_curr(:,:,j_marker,j_state);
                                obj.copycurrtoprops();
                                obj.rejected(j_up) = obj.rejected(j_up) + 1/(obj.nomarkers*obj.nostates); % one more rejected proposal
                            end
                            obj.samplecounter(j_up) = obj.samplecounter(j_up) + 1/(obj.nomarkers*obj.nostates);
                        end     % end of markers loop
                    end     % end of states loop
                elseif( strcmp(obj.errorupdatetype,'Gibbs') )
                    % get new samples for all  error terms:
                    obj.updateerrorprop_Gibbs();         % gives obj.merr_prop
                    obj.logtarget_prop = obj.getabsolutetargetprobability(-1,-1);
                    % accept directly:
                    obj.copypropstocurr();
                else    % i.e. unknown sampler
                    fprintf( ' Warning - updatecurrentsample: Unkown errorupdatetype of %s\n', obj.errorupdatetype );
                    return;
                end     % end of going through different kinds of samplers for error

                % check new perspective proposal:
                % ...get new update for perspectives sequentially for each datapoint:
                j_up = j_up + 1;
                if( obj.without>=2 )
                    fprintf( ' Info - updatecurrentsample (%s)(%d): Perspectives\n', obj.comment,obj.MCit );
                    obj.comparecurrentandproposedpars();
                end
                for j_data = 1:obj.nodata                       % sequential update for each sister-pair individually
                    if( ~all(~obj.posstatevskin(1:obj.nostates,j_data)) )	% i.e. if not all states impossible
                        % ....get new proposal:
                        obj.updateperspectiveprops(j_data, j_up);
                        obj.updatepositionpredictions(j_data,-1);	% not yet done in updateperspecitveprops.
                        [logtargetratio,obj.logtarget_prop] = obj.gettargetprobabilityratio(j_data,-1);
                        if( rand()<exp(logtargetratio) )      	% accept
                            if( (obj.without>=3) && (j_data==1) )
                                fprintf( '%d: accept persp_prop at j_data=%3d, (%1.3e)(%1.6e-->%1.6e)\n', obj.MCit, j_data, logtargetratio,obj.logtarget_curr,obj.logtarget_prop );
                            end
                            % ....copy proposals into current samples:
                            obj.copypropstocurr();
                        else                                	% reject
                            if( (obj.without>=3) && (j_data==1) )
                                fprintf( '%d: reject persp_prop at j_data=%3d, (%1.3e)(%1.6e-->%1.6e)\n', obj.MCit, j_data, logtargetratio,obj.logtarget_curr,obj.logtarget_prop );
                            end
                            obj.copycurrtoprops();
                            obj.rejected(j_up) = obj.rejected(j_up) + 1/(obj.nodata);   % one more rejected proposal
                        end     % end of accept/reject step
                        obj.samplecounter(j_up) = obj.samplecounter(j_up) + 1/(obj.nodata);
                    end     % end if some state possible for this data
                end     % end of sequential update of perspectives
                % update translations as Gibbs proposal:
                if( 1 )
                    if( obj.without>=2 )
                        fprintf( ' Info - updatecurrentsample (%s)(%d): Translations\n', obj.comment,obj.MCit );
                        obj.comparecurrentandproposedpars();
                    end
                    obj.updatetranslations_Gibbs(-1);
                    obj.updatepositionpredictions(-1,-1);
                    obj.logtarget_prop = obj.getabsolutetargetprobability(-1,-1);
                    % ...accept for sure:
                    obj.copypropstocurr();
                end     % end if 1
                    
                % state update:
                if( obj.nostates>1 )
                    % ...state fractions:
                    if( obj.without>=2 )
                        fprintf( ' Info - updatecurrentsample (%s)(%d): Statefrac\n', obj.comment,obj.MCit );
                        obj.comparecurrentandproposedpars();
                    end
                    obj.updatestatefrac_Gibbs();            % get Gibbs proposal
                    obj.logtarget_prop = obj.getabsolutetargetprobability(-1,-1);
                    % accept for sure:
                    obj.copypropstocurr();
                    % ...occupation numbers:
                    if( obj.without>=2 )
                        fprintf( ' Info - updatecurrentsample (%s)(%d): Stateocc\n', obj.comment,obj.MCit );
                        obj.comparecurrentandproposedpars();
                    end
                    j_up = j_up + 1;
                    order = randperm(obj.nodata);           % random order
                    stateoptions = 1:obj.nostates;          % all states possible
                    for j_data = order
                        stateoptions_here = stateoptions( obj.posstatevskin(:,j_data) );    % those states allowed here
                        nostateoptions_here = numel(stateoptions_here);
                        if( nostateoptions_here>1 )         % i.e. more than one state allowed for this datapoint
                            newstate = ~obj.statevskin_prop(stateoptions_here,j_data); newstate = stateoptions_here(newstate); newstate = newstate(ceil(rand()*(nostateoptions_here-1)));
                            obj.statevskin_prop(:,j_data) = false(obj.nostates,1);  obj.statevskin_prop(newstate,j_data) = (1==1);  % random new state, but different from current one
                            obj.keepmeansinvariant( j_data );	% correct to keep mean invariant
                            % ...get targetdensity:
                            obj.updatepositionpredictions(j_data,-1);
                            [logtargetratio,obj.logtarget_prop] = obj.gettargetprobabilityratio(j_data,-1);
                            % ...decide to accept/reject:
                            if( rand()<exp(logtargetratio) )	% accept
                                if( obj.without>=3 )
                                    fprintf( ' Info - updatecurrentsample (%d)(%s): accept state of j_data=%d (%d-->%d). (logtargetratio = %1.3e)\n', obj.MCit,obj.comment, j_data, stateoptions(obj.statevskin_curr(:,j_data)),stateoptions(obj.statevskin_prop(:,j_data)), logtargetratio );
                                end
                                obj.copypropstocurr();
                            else                            % reject
                                if( obj.without>=3 )
                                    fprintf( ' Info - updatecurrentsample (%d)(%s): reject state of j_data=%d (%d-->%d). (logtargetratio = %1.3e)\n', obj.MCit,obj.comment, j_data, stateoptions(obj.statevskin_curr(:,j_data)),stateoptions(obj.statevskin_prop(:,j_data)), logtargetratio );
                                end
                                obj.copycurrtoprops();
                                obj.rejected(j_up) = obj.rejected(j_up) + 1/obj.nodata;
                            end     % end of accept reject
                            obj.samplecounter(j_up) = obj.samplecounter(j_up) + 1/obj.nodata;
                        end     % end if multiple different states allowed for this datapoint
                    end     % end of going through all datapoints
                end     % end if multiple states
            end     % end of subsampling
            
            %obj.copycurrtoprops();
            % update history:
            obj.flpos_hist(:,:,:,obj.MCit) = obj.flpos_curr;
            obj.merr_hist(:,:,:,:,obj.MCit) = obj.merr_curr;
            obj.trans_hist(:,:,obj.MCit) = obj.trans_curr;
            obj.rot_hist(:,:,:,obj.MCit) = obj.rot_curr;
            obj.statefrac_hist(:,:,obj.MCit) = obj.statefrac_curr;
            obj.statevskin_hist(:,:,obj.MCit) = obj.statevskin_curr;
            obj.logtarget_hist(1,obj.MCit) = obj.logtarget_curr;
            obj.writestatestext();  % write to external file
        end     % end of updatecurrentsample function
        
        function[] = updatetruegeometryprops( obj, j_marker, j_state,j_up )
            % updates the template positions via random walk Metropolis Hastings
            
            % set auxiliary parameters:
            staterange = j_state;
            if( j_state<0 )
                staterange = 1:obj.nostates;
            end
            markerrange = j_marker;
            if( j_marker<0 )
                markerrange = 1:obj.nomarkers;
            end
                    
            for j_state = staterange
                for j_marker = markerrange
                    obj.flpos_prop(:,j_marker,j_state) = obj.flpos_curr(:,j_marker,j_state) + obj.stepsizes(j_up)*obj.sampleisoball();
                    if( j_marker==1 )   % first fluorophore special, because centre of rotation
                        obj.updatepositionpredictions( -1,j_state );
                        obj.trans_prop = obj.trans_curr - squeeze(obj.flpos_pred_prop(:,2,:)-obj.flpos_pred_curr(:,2,:));   % keep other fluorphores where they were
                        obj.updatepositionpredictions( -1,j_state );
                    end     % end if first fluorophore 
                end     % end of markerrange loop
            end     % end of staterange loop
        end     % end of updatetruegeometryprops function
        
        function[] = updateerrorprop( obj, j_marker, j_state, j_up )
            % updates the measurement error standard deviation merr1, merr2 or merr3
            % via rwMH update with Gaussian restricted to positive values
            
            % set auxiliary parameters:
            staterange = j_state;
            if( j_state<0 )
                staterange = 1:obj.nostates;
            end
            markerrange = j_marker;
            if( j_marker<0 )
                markerrange = 1:obj.nomarkers;
            end
            
            for j_state = staterange
                for j_marker = markerrange
                    obj.merr_prop(1,1,j_marker,j_state) = abs( obj.merr_curr(1,1,j_marker,j_state) + obj.stepsizes(j_up)*(2*rand()-1) );
                    obj.merr_prop(2,2,j_marker,j_state) = obj.merr_prop(1,1,j_marker,j_state);  % x=y
                    obj.merr_prop(3,3,j_marker,j_state) = abs( obj.merr_curr(3,3,j_marker,j_state) + obj.stepsizes(j_up)*(2*rand()-1) );
                end     % end of markerrange loop
            end     % end of staterange loop
        end     % end of updateerrorprop function
        function[] = updateerrorprop_Gibbs( obj )
            % updates the measurement error merr via Gibbs sampler
            
            for j_state = 1:obj.nostates
                % get applicable datapoints:
                j_data = obj.statevskin_curr(j_state,:);
                ntotal = sum(j_data);
                for j_marker = 1:obj.nomarkers  % go through each fluorophore separately
                    % priors:
                    sp = obj.shapepar_prior{j_state}(1);                          	% '1-3/2' is flat on sigma; '-3/2' for flat prior on sigma; '1' for getting rid of 1 in gamma-distr. definition; 
                    rp = (2/sp)./(([2,0,2].*obj.merrmeans_prior{j_state}(j_marker,:)).^2)';% '[2,0,2].*' to get mean of std right, not precision
                    % get parameters in gamma-distribution:
                    shapeparxy = sp + 0.5*2*ntotal;                                 % '0.5' to get sqrt(tau) in front; '2' for two space dimensions
                    shapeparz = sp + 0.5*1*ntotal;                                  % '0.5' to get sqrt(tau) in front; '1' for one space dimensions
                    ratepar = 1./rp + zeros(3,1);                                   % initialise according to priors
                    ratepar = ratepar + sum( (obj.flpos_data(:,j_marker,j_data)-obj.flpos_pred_curr(:,j_marker,j_data)).^2, 3);  	% fl1; sum over all datapoints of that state
                    rateparxy = ratepar(1)+ratepar(2);  rateparz = ratepar(3);      % sum over space dimensions
                    keepontrying = 0;                                               % to avoid singular results
                    while( keepontrying>=0 )
                        obj.merr_prop(1,1,j_marker,j_state) = gamrnd(shapeparxy,1/(rateparxy/2));% sample precision = 1/sigma^2
                        obj.merr_prop(1,1,j_marker,j_state) = 1/sqrt(obj.merr_prop(1,1,j_marker,j_state));  % get actual standard deviation
                        obj.merr_prop(2,2,j_marker,j_state) = obj.merr_prop(1,1,j_marker,j_state);  % x=y
                        obj.merr_prop(3,3,j_marker,j_state) = gamrnd(shapeparz,1/(rateparz/2));	% sample precision = 1/sigma^2
                        obj.merr_prop(3,3,j_marker,j_state) = 1/sqrt(obj.merr_prop(3,3,j_marker,j_state));  % get actual standard deviation
                        if( all(~isinf(obj.merr_prop(:,:,j_marker,j_state)),'all') && all(obj.merr_prop(:,:,j_marker,j_state)>=0,'all') )
                            keepontrying = -1;
                        else
                            keepontrying = keepontrying + 1;
                            if( keepontrying>100 )
                                fprintf( ' Warning - updateerrorprop_Gibbs (%d)(%s): Cannot find non-singular error for st=%d, fl=%d:\n', obj.MCit,obj.comment, j_state,j_marker );
                                disp(shapeparxy), disp(shapeparz), disp(ratepar), disp(obj.merr_prop);
                                return;
                            end     % end if tried too long
                        end     % end if pathological
                    end     % end of keeptrying to get non-singular result
                end     % end of markers loop
            end     % end of state loop
        end     % end of updateerrorprop_Gibbs function
        
        function[] = updateperspectiveprops( obj, j_data, j_up )
            % updates perspective proposals for each individual sister-sister-pair with a rwMH step on the real space and then transforms them into Euler-angles
            
            if( j_data<0 )       % go through all data indices
                for j_data = 1:obj.nodata
                    obj.updateperspectiveprops(j_data, j_up);
                end     % end of going through data
         	else                    % only update the data index given
                stepsizehere = 1/sqrt(sum(obj.stepsizes((obj.nostates-1)*obj.nomarkers +(1:obj.nomarkers)).^(-2)));   % all markers of this state combined
                obj.trans_prop(:,j_data) = obj.trans_curr(:,j_data) + stepsizehere*obj.sampleisoball(); % same steps as for fluorophores
                isotropicrandomref = obj.sampleisosphere([],2);
                perturbedref = obj.sampleisosphere(isotropicrandomref,1-cos(min(pi,obj.stepsizes(j_up))));
                a = cross(isotropicrandomref,perturbedref); alpha = acos(isotropicrandomref'*perturbedref);
                obj.rot_prop(:,:,j_data) = obj.get3drotation( a, alpha )*obj.rot_curr(:,:,j_data);
            end     % end distinguishing dataindex
        end     % end of updateperspectiveprops4 function
        function[] = updatetranslations_Gibbs( obj, j_data )
            % updates translations of each triangle sequentially via Gibbs sampler
            % scalefactor is to use larger independent proposal
            
            % get auxiliary parameters:
            stateoptions = 1:obj.nostates;  % all states possible
            tau = zeros(3,3,obj.nomarkers,obj.nostates); sotau = zeros(3,3,obj.nostates);    isotau = zeros(3,3,obj.nostates);
            for j_state = 1:obj.nostates
                for j_marker = 1:obj.nomarkers
                    tau(:,:,j_marker,j_state) = diag( diag(obj.merr_curr(:,:,j_marker,j_state)).^(-2) );    % uses diagonal structure
                    sotau(:,:,j_state) = sotau(:,:,j_state) + tau(:,:,j_marker,j_state);
                end     % end of markers loop
                isotau(:,:,j_state) = sotau(:,:,j_state)^(-1);    % (inverse of) sum of precision matrices
            end     % end of states loop
            if( j_data==-1 )
                datarange = 1:obj.nodata; datarange = datarange(sum(obj.posstatevskin,1)>0);  % i.e. only those states that do contribute to at least one datapoint
            else
                datarange = j_data;
            end     % end if dataindex is given
            
            % go through each datapoint:
            for j_data = datarange
                j_state = stateoptions( obj.statevskin_curr(:,j_data) );
                mus = zeros(3,1);       % initialise means for this datapoint
                for j_marker = 1:obj.nomarkers
                    % ...get predictions without translation:
                    flpos_rot = obj.rot_prop(:,:,j_data)*(obj.flpos_curr(:,j_marker,j_state) - obj.flpos_curr(:,1,j_state) ) + obj.flpos_curr(:,1,j_state);
                    % compare to observations, to get mean:
                    mus = mus + ( tau(:,:,j_marker,j_state)*(obj.flpos_data(:,j_marker,j_data) - flpos_rot) );
                end     % end of markers loop
                mus = isotau(:,:,j_state)*mus;
                % sample from normal distribution:
                obj.trans_prop(:,j_data) = mvnrnd( mus', isotau(:,:,j_state)  )';
            end     % end of going through each datapoint
        end     % end of updatetranslations_Gibbs function
        
        function[ rejectdirectly ] = updatestatefrac_rwMH( obj, j_state,j_dset, j_up )
            % updates statefrac_prop
            % random walk proposal for j_state, by picking a random other state to keep the sum constant
            
            % set auxiliary parameters:
            rejectdirectly = 0;         % don't reject directly; only if beyond [0,1]-range
            jj_state = 1:obj.nostates; jj_state(j_state) = []; jj_state = jj_state(ceil((obj.nostates-1)*rand()));	% choose state among the remaining ones at random
            
            % get new proposal for j_state:
            change = obj.stepsizes(j_up)*(2*rand()-1);
            obj.statefrac_prop(j_state,j_dset) = obj.statefrac_curr(j_state,j_dset) + change;
            % keep sum constant:
            obj.statefrac_prop(jj_state,j_dset) = obj.statefrac_curr(jj_state,j_dset) - change;   % sum of jj_state and j_state must remain constant
            if( (min(obj.statefrac_prop)<0) || (max(obj.statefrac_prop)>1) )
                rejectdirectly = 1;
            end     % end if rejectdirectly
        end     % end of updatestatefrac_rwMH function
        function[ loghastingsterm ] = updatestatefrac_Gibbs( obj )
            % updates the state fraction based on Dirichlet distribution
            % not a Gibbs-sampler if more than some datapoint are not allowed
            % to have all states but more than one (i.e is Gibbs, if not
            % more than two states or all datapoints allowed all states or
            % only one state)
            
            % get auxiliary parameters:
            loghastingsterm = 0;    % initialise
            for j_dset = 1:obj.nodsets
                morethanonestatepossible = (sum(obj.posstatevskin,1)>1)&(obj.datasetvskin(j_dset,:));% logical '1' if more than one state possible for datapoint at respective index
                if( isempty(morethanonestatepossible) )
                    nstate = zeros(obj.nostates,1); % indifferent, if no datapoints allow multiple states
                else
                    nstate = sum( obj.statevskin_curr(:,morethanonestatepossible), 2 );    % column vector with the occupation numbers per state
                end     % end if no datapoints with multiple states
                
                % first generate Gamma-distributed random numbers:
                keepondoing = (1==1);
                while( keepondoing )
                    buffer = gamrnd( nstate+1, 1 );                 	% intermediate buffer, that still needs to be normalised
                    if( min(buffer)>0 ), keepondoing = (0==1);   end    % found non-zero state fractions
                end     % end of keeping searching for finite probabilities
                
                % transform into Dirichlet random varialbes:
                obj.statefrac_prop(:,j_dset) = buffer/sum(buffer);

                % get loghastingsterm:
                loghastingsterm = loghastingsterm + (log(obj.statefrac_curr(nstate>0,j_dset))-log(obj.statefrac_prop(nstate>0,j_dset)))'*nstate(nstate>0);  % prop of curr/prop of prop
            end     % end of datasets loop
        end     % end of updatestatefrac_Gibbs function
        
        function[ logtargetratio,logtargetprop ] = gettargetprobabilityratio( obj, j_data, j_state, withpriors )
            % computes log of targetratio
   
            % set auxiliary parameters:
            logtargetratio = 0;                     % initialise log of targetratio
            dataoptions = 1:obj.nodata;             % all datapoints possible
            stateoptions = 1:obj.nostates;          % all states possible
            if( (j_data<0) && (j_state<0) )         % both are full range
                datarange = dataoptions;
                staterange = stateoptions;
            elseif( j_data<0 )                      % state specified but data full range
                datarange = [dataoptions( obj.statevskin_curr(j_state,:) ),dataoptions( obj.statevskin_prop(j_state,:) )];
                datarange = unique(datarange);
                staterange = j_state;
            elseif( j_state<0 )                     % data specified but state full range
                datarange = j_data;
                staterange = [stateoptions( obj.statevskin_curr(:,j_data) ),stateoptions( obj.statevskin_prop(:,j_data) )];
                staterange = unique(staterange);
            else                                    % both specified
                datarange = j_data;
                staterange = j_state;
            end     % end if range
            t_prop_l = zeros(3,3,obj.nomarkers,obj.nostates);   t_curr_l = zeros(3,3,obj.nomarkers,obj.nostates);
            for j_state = staterange
                for j_marker = 1:obj.nomarkers
                    t_prop_l(:,:,j_marker,j_state) = diag(diag( obj.merr_prop(:,:,j_marker,j_state).^(-1)) );    % for easier computation (1/stddev); uses diagonal shape
                    t_curr_l(:,:,j_marker,j_state) = diag(diag( obj.merr_curr(:,:,j_marker,j_state).^(-1)) );    % for easier computation (1/stddev); uses diagonal shape
                end     % end of markers loop
            end     % end of going through states
            if( nargin<4 )  % i.e. withpriors not given
                withpriors = 1; % i.e. include priors in computation
            end     % end of filling in withpriors, if necessary
            
            % check priors:
            logpriorratio = 0;           % initialise log of prior
           	if( withpriors==1 )
                % check hard boundaries first: (onyl on prop)
                % ...for errors:
                for j_state = staterange
                    if( ~all(~isinf(obj.merr_prop(:,:,:,j_state)),'all')  )
                        logtargetratio = -Inf;  logtargetprop = -Inf; return;
                    end     % end if errcond
                end
                % ...for occupation numbers:
                nstates_prop = sum(obj.statevskin_prop,2);  % total occupation numbers
                if( min(nstates_prop)<3 )                   % want to have at least three in each state
                  	logtargetratio = -Inf;  logtargetprop = -Inf; return;
                end     % end if state completely unoccupied
                % take length prior into account:
                logpriorratio = logpriorratio + obj.getabsoluteloglengthprior( obj.flpos_prop )-obj.getabsoluteloglengthprior( obj.flpos_curr );
                % take prior on t_prop,t_curr into account:
                for j_state = staterange        % others are zero anyways
                    if( (min(obj.merr_prop(:,:,:,j_state),[],'all')<0) )
                        logtargetratio = -Inf;  logtargetprop = -Inf; return;
                    end     % end if negative measurement errors
                    for j_marker = 1:obj.nomarkers
                        merrmeansh2prior = diag(obj.merrmeans_prior{j_state}(j_marker,:).^2)*(obj.shapepar_prior{j_state}(j_marker));
                        logpriorratio = logpriorratio + log(det(t_prop_l(:,:,j_marker,j_state)^2))*( obj.shapepar_prior{j_state}(1)-(1-(3/2)) ) - trace( merrmeansh2prior*(t_prop_l(:,:,j_marker,j_state)^2) );   % note: this is the relative prior to the flat Lebesgue measure on merr1 (i.e. Radon-Nikodym derivative wrt to this), which has shapepar (1-(3/2)) and ratepar [0,0,0] for the three dimensions
                        logpriorratio = logpriorratio - log(det(t_curr_l(:,:,j_marker,j_state)^2))*( obj.shapepar_prior{j_state}(1)-(1-(3/2)) ) - trace( merrmeansh2prior*(t_curr_l(:,:,j_marker,j_state)^2) );
                    end     % end of markers loop
                end     % end of state loop
            end     % end if withpriors
            
            % check likelihood:
            dsets = 1:obj.nodsets;
            for j_data = datarange
                j_state_prop = stateoptions(obj.statevskin_prop(:,j_data));
                j_state_curr = stateoptions(obj.statevskin_curr(:,j_data));
                j_dset = dsets(obj.datasetvskin(:,j_data));
                for j_marker = 1:obj.nomarkers
                    % get state and corresponding measurement error (right sister):
                    t_prop = t_prop_l(:,:,j_marker,j_state_prop);      % for easier computation (1/stddev)
                    t_curr = t_curr_l(:,:,j_marker,j_state_curr);      % for easier computation (1/stddev)
                    % get target-value (right sister):
                    fct = (t_prop)*( obj.flpos_pred_prop(:,j_marker,j_data)-obj.flpos_data(:,j_marker,j_data) );
                    logtargetratio = logtargetratio + (-1/2)*sum(fct.^2) + log(det(t_prop));
                    fct = (t_curr)*( obj.flpos_pred_curr(:,j_marker,j_data)-obj.flpos_data(:,j_marker,j_data) );
                    logtargetratio = logtargetratio - (-1/2)*sum(fct.^2) - log(det(t_curr));
                end     % end of markers loop
                if( obj.nostates>1 )
                    logtargetratio = logtargetratio + log(obj.statefrac_prop(j_state_prop,j_dset)/sum(obj.statefrac_prop(obj.posstatevskin(:,j_data),j_dset)));
                    logtargetratio = logtargetratio - log(obj.statefrac_curr(j_state_curr,j_dset)/sum(obj.statefrac_curr(obj.posstatevskin(:,j_data),j_dset)));
                end     % end if multistate label model
            end     % end of data loop
            % take prior into account:
            logtargetratio = logtargetratio + logpriorratio;
            logtargetprop = obj.logtarget_curr + logtargetratio;
            % check if pathological:
            if( isnan(logtargetratio) )
                fprintf( ' Info - gettargetratio (%d)(%s): Num. probl. for logtargetratio = %1.3e\n', obj.MCit,obj.comment, logtargetratio );
                logtargetratio = -Inf; logtargetprop = -Inf;
            end     % end if pathological
        end     % end of gettargetratio function
        function[ logtarget ] = getabsolutetargetprobability( obj, j_data, j_state, withpriors )
            % computes (for proposed parameters) targetdensity (up to a constant) wrt to Lebesgue measure of positions, stddevs, translations and uniform measure on sphere
            %fprintf( ' Info - getabsolutetarget (%d)(%s):\n', obj.MCit,obj.comment ); 
            
            % set auxiliary parameters:
            logtarget = 0;                          % initialise log of target density
            dataoptions = 1:obj.nodata;             % all datapoints possible
            stateoptions = 1:obj.nostates;          % all states possible
            if( (j_data<0) && (j_state<0) )         % both are full range
                datarange = dataoptions;
                staterange = stateoptions;
            elseif( j_data<0 )                      % state specified but data full range
                datarange = dataoptions( obj.statevskin_prop(j_state,:) );
                staterange = j_state;
            elseif( j_state<0 )                     % data specified but state full range
                datarange = j_data;
                staterange = stateoptions( obj.statevskin_prop(:,j_data) );
            else                                    % both specified
                datarange = j_data;
                staterange = j_state;
            end     % end if range
            t_prop_l = zeros(3,3,obj.nomarkers,obj.nostates);
            for j_state = staterange
                for j_marker = 1:obj.nomarkers
                    t_prop_l(:,:,j_marker,j_state) = diag(diag( obj.merr_prop(:,:,j_marker,j_state).^(-1)) );    % for easier computation (1/stddev); uses diagonal shape
                end     % end of markers loop
            end     % end of going through states
            if( nargin<4 )  % i.e. withpriors not given
                withpriors = 1; % i.e. include priors in computation
            end     % end of filling in withpriors, if necessary
            
            % check priors:
            logprior = 0;           % initialise log of prior
           	if( withpriors==1 )
                % check hard boundaries first:
                % ...for errors:
                for j_state = staterange
                    if( ~all(~isinf(obj.merr_prop(:,:,:,j_state)),'all')  )
                        logtarget = -Inf; return;
                    end     % end if errcond
                end
                % ...for occupation numbers:
                nstates_prop = sum(obj.statevskin_prop,2);  % total occupation numbers
                if( min(nstates_prop)<3 )                   % want to have at least three in each state
                  	logtarget = -Inf; return;
                end     % end if state completely unoccupied
                % take length prior into account:
                logprior = logprior + obj.getabsoluteloglengthprior( obj.flpos_prop );
                % take prior on t_prop into account:
                for j_state = staterange        % others are zero anyways
                    if( (min(obj.merr_prop(:,:,:,j_state),[],'all')<0) )
                        logtarget = -Inf; return;
                    end     % end if negative measurement errors
                    for j_marker = 1:obj.nomarkers
                        merrmeansh2prior = diag(obj.merrmeans_prior{j_state}(j_marker,:).^2)*(obj.shapepar_prior{j_state}(j_marker));
                        logprior = logprior + log(det(t_prop_l(:,:,j_marker,j_state)^2))*( obj.shapepar_prior{j_state}(1)-(1-(3/2)) ) - trace( merrmeansh2prior*(t_prop_l(:,:,j_marker,j_state)^2) );   % note: this is the relative prior to the flat Lebesgue measure on merr1 (i.e. Radon-Nikodym derivative wrt to this), which has shapepar (1-(3/2)) and ratepar [0,0,0] for the three dimensions 
                    end     % end of markers loop
                end     % end of state loop
            end     % end if withpriors
            
            % check likelihood:
            dsets = 1:obj.nodsets;
            for j_data = datarange
                j_state_prop = stateoptions(obj.statevskin_prop(:,j_data));
                j_dset = dsets(obj.datasetvskin(:,j_data));
                for j_marker = 1:obj.nomarkers
                    % get state and corresponding measurement error (right sister):
                    t_prop = t_prop_l(:,:,j_marker,j_state_prop);      % for easier computation (1/stddev)
                    % get target-value (right sister):
                    fct = (t_prop)*( obj.flpos_pred_prop(:,j_marker,j_data)-obj.flpos_data(:,j_marker,j_data) );
                    logtarget = logtarget + (-1/2)*sum(fct.^2) + log(det(t_prop));
                end     % end of markers loop
                if( obj.nostates>1 )
                    logtarget = logtarget + log(obj.statefrac_prop(j_state_prop,j_dset)/sum(obj.statefrac_prop(obj.posstatevskin(:,j_data),j_dset)));
                end     % end if multistate label model
            end     % end of data loop
            % take prior into account:
            logtarget = logtarget + logprior;
            % check if pathological:
            if( isnan(logtarget) )
                fprintf( ' Info - getabsolutetarget (%d)(%s): Num. probl. for logtarget = %1.3e\n', obj.MCit,obj.comment, logtarget );
                logtarget = -Inf;
            end     % end if pathological
        end     % end of getabsolutetarget function
        function[ loglengthprior ] = getabsoluteloglengthprior( obj, flpos,loglengthsmatrix )
            % outputs the prior of the proposed positions
            
            loglengthprior = 0;
            if( obj.lengthpriortype==0 )    % flat prior on true positions
            	% keep as is; constant prior on true positions (not lengths)
            elseif( obj.lengthpriortype==1 )% length dependent
                if( (nargin<3) || isempty(loglengthsmatrix) )
                    % get lengths:
                    nolengths = (obj.nomarkers-1)*obj.nomarkers/2;
                    lengthlist = zeros(nolengths,obj.nostates);     % only needed at top level, i.e. when lengthsmatrix is not established, yet
                    loglengthsmatrix{obj.nostates} = zeros(obj.nomarkers,obj.nomarkers);
                    for j_state = 1:obj.nostates
                        j_length = 0;
                        loglengthsmatrix{j_state} = zeros(obj.nomarkers,obj.nomarkers);   % initialise
                        for j_marker = 1:(obj.nomarkers-1)
                            for jj_marker = (1+j_marker):obj.nomarkers
                                j_length = j_length + 1;
                                lengthlist(j_length,j_state) = sqrt(sum((flpos(:,j_marker,j_state)-flpos(:,jj_marker,j_state)).^2,1));
                                loglengthsmatrix{j_state}(j_marker,jj_marker) = log(lengthlist(j_length,j_state)); % don't need second half
                            end     % end of inner marker loop
                        end     % end of outer markers loop
                    end     % end of states loop
                    % apply given length prior (relative to Lebesgue measure):
                    if( obj.lengthpriordensity.type==1 )    % intervals
                        for j_state = 1:obj.nostates
                            cond = all( abs(obj.lengthpriordensity.pars{j_state}(:,1)-lengthlist(:,j_state))<=obj.lengthpriordensity.pars{j_state}(:,2) );
                            if( ~cond )         % condition violated/ proposals outside of prior intervals
                                loglengthprior = -Inf; return;  % outside of intervals; don't need to check other states
                            end     % end if violates hard boundary
                        end     % end of state loop
                    end     % end of checking density
                end     % end of getting lengthsmatrix myself
                % get minimal legnth:
                minloglength = +Inf(3,obj.nostates);        % [loglength,j_marker,jj_marker]x[j_state]
                nomarkershere = size(loglengthsmatrix{1},1);% might be shrunk in recursive relation
                for j_state = 1:obj.nostates
                    for j_marker = 1:(nomarkershere-1)
                        for jj_marker = (j_marker+1):nomarkershere
                            if( loglengthsmatrix{j_state}(j_marker,jj_marker)<minloglength(1,j_state) )
                                minloglength(:,j_state) = [ loglengthsmatrix{j_state}(j_marker,jj_marker); j_marker; jj_marker ];
                            end     % end of found new minimal length
                        end     % end of inner markers loop
                    end     % end of outer markers loop
                end     % end of states loop
                % get prior recursively:
                loglengthprior = loglengthprior - 2*sum(minloglength(1,:),2);
                if( nomarkershere>2 )   % more than two markers left, otherwise at bottom of recursion
                    loglengthsmatrix_sub1 = loglengthsmatrix;   loglengthsmatrix_sub2 = loglengthsmatrix;
                    for j_state = 1:obj.nostates
                        loglengthsmatrix_sub1{j_state}(minloglength(2,j_state),:) = []; loglengthsmatrix_sub1{j_state}(:,minloglength(2,j_state)) = [];  % delte datapoints, that are endpoints of current minlength
                        loglengthsmatrix_sub2{j_state}(minloglength(3,j_state),:) = []; loglengthsmatrix_sub2{j_state}(:,minloglength(3,j_state)) = [];
                    end     % end of states loop
                    loglengthprior = loglengthprior + ( obj.getabsoluteloglengthprior([],loglengthsmatrix_sub1)+obj.getabsoluteloglengthprior([],loglengthsmatrix_sub2) )/2;    % no positions needed, once lengths are calculated
                end     % end if down to two markers
            else    % unknown type
                fprintf( ' (%s) Warning - getabsoluteloglengthprior (%d): Unknown lengthpriortype %d.\n', obj.comment,obj.MCit, obj.lengthpriortype );
            end     % end of distinguishing lengthpriortypes
        end     % end of getabsoluteloglengthprior function
        
        function[] = updatepositionpredictions( obj, j_data,j_state )
            % computes the predicted, proposed positions of the fluorophores
            
            % set auxiliary parameters:
            if( j_state<0 )
                if( obj.nostates>1 )
                    kinsofthesestates = (~all(~obj.statevskin_prop(:,:),1));   % should be all ones; all acts column-wise
                else
                    kinsofthesestates = obj.posstatevskin(1,:);
                end     % end if multiple states in mixture
            else
                kinsofthesestates = obj.statevskin_prop(j_state,:);
            end
            stateoptions = 1:obj.nostates;  % all states possible
            
            if( j_data<0 )       % go through all data
                datarange = 1:obj.nodata; datarange = datarange(kinsofthesestates);                             % only those datapoints that have given state
                for j_data = datarange
                    obj.updatepositionpredictions(j_data,j_state);
                end     % end of data loop
            else                    % only update the data index given
                if( kinsofthesestates(j_data)==1 )      % otherwise don't update anything
                    j_state = stateoptions(obj.statevskin_prop(:,j_data));
                    for j_marker = 1:obj.nomarkers
                        obj.flpos_pred_prop(:,j_marker,j_data) = obj.rot_prop(:,:,j_data)*(obj.flpos_prop(:,j_marker,j_state)-obj.flpos_prop(:,1,j_state)) + obj.flpos_prop(:,1,j_state) + obj.trans_prop(:,j_data);
                    end     % end of markers loop
                end     % end of checking compatibility of stateindex and dataindex
            end     % end of distinguishing datapoints
        end     % end of updatepositionpredictions function
        function[ R3d ] = get3drotation( ~, a,alpha )
            % gets rotation matrix that rotates the angle phi around the given axis a
            
            na = norm(a);
            if( na>0 )              	% distinction avoids proplems with a not being normalised
                % normalise the axis:
                a = a(:)/na;            % normalised column vector
                % get 3d rotation matrix around a:
                R3d = eye(3)*cos(alpha) + a*(a')*(1-cos(alpha)) + [[0,-a(3),a(2)];[a(3),0,-a(1)];[-a(2),a(1),0]]*sin(alpha);
                %{
                if( sum(abs(imag(R3d(:))))>1e-20 )
                    fprintf( ' Warning - get3drotation: Large imaginary parts in R3d: %1.3e\n', sum(sum(imag(R3d).^2)) );
                    fprintf( ' %+1.3e, %+1.3e, %+1.3e\n', R3d(1,:) );
                    fprintf( ' %+1.3e, %+1.3e, %+1.3e\n', R3d(2,:) );
                    fprintf( ' %+1.3e, %+1.3e, %+1.3e\n', R3d(3,:) );
                end
                %}
                R3d = real(R3d);        % round off wrong imaginary terms
            else
                R3d = eye(3);       	% identity; likely a==0, as crossproduct for alpha=0
            end
        end     % end of get3drotation function
        function[ isosphrnd ] = sampleisosphere( obj, meanaxis, perturbation )
            % samples isotropically from surface of unit sphere in 3d
            % meanaxis is symmetry axis (must be on unitsphere); can be empty if on entire sphere
            % perturbation is in [0,2] and measured along the meanaxis (i.e. 1-cos(phi))
            % outputs a columnvector
            
            if( isempty(meanaxis) )     % just uniform on entire sphere (should be slightly faster than more general case)
                xy = 2*pi*rand(); z = 2*rand()-1; nz = sqrt(1-z^2);
                isosphrnd = [ nz*cos(xy); nz*sin(xy); z ];
            else                        % local perturbation around meanaxis
                xy = 2*pi*rand(); z = 1-perturbation*rand(); nz = sqrt(1-z^2);
                r = [nz*cos(xy);nz*sin(xy);z];
                a = cross([0;0;1],meanaxis); alpha = acos([0;0;1]'*meanaxis);
                isosphrnd = obj.get3drotation(a,alpha)*r;
            end
        end     % end of sampleisosphere function
        function[ isoballrnd ] = sampleisoball( obj )
            % samples from a uniform unitball in 3d
            % output is columnvector
            
            isoballrnd = (rand()^(1/3)) * obj.sampleisosphere([],2);
        end     % end of sampleisoball function
        function[ R3d ] = sampleisorotation( obj )
            isosphrnd = obj.sampleisosphere([],2);	% isotropic on sphere
            keeptrying = (1==1);
            while( keeptrying )
                angle = 2*pi*rand();
                if( rand()<((1-cos(angle))/2) )     % rejection sampler for rotation angle
                    keeptrying = (0==1);            % succeded
                end     % end of rejection sampler for getting rotation angle
            end     % end of keeptrying
            R3d = obj.get3drotation( isosphrnd,angle );
        end     % end of sampleisorotation function
        function[] = optimisestepsizes( obj )
            % adjust stepsizes during burnin, so rejectionrates are within reasonable ranges
            % to be called after parameter update of current MCit
            
            % auxiliary parameters:
            adjfactor = 2;                          % should be >=1
            middle = 0.775; span = 0.075;           % 'reasonable' range for rejection rates
            optstepsfreq = max( 10, ceil(sqrt(obj.burnin)) );	% how often optimisestepsizes becomes active
            minsamples = 250;                      	% minimum number of samples to adjust stepsize
            
            if( obj.MCit==1 )
                obj.nextoptstepstime = obj.burnin - floor(obj.burnin/optstepsfreq)*optstepsfreq + optstepsfreq;
            end     % end of setting first optstepstime
            if( (obj.MCit==obj.nextoptstepstime) && (obj.MCit<=obj.burnin) )	% only optimise during burnin phase
                % check each update separately:
                for j_up = 1:numel(obj.stepsizes)
                    if( obj.samplecounter(j_up)>=minsamples )
                        rejectionratesincelastadjustment = obj.rejected(j_up)/obj.samplecounter(j_up);
                        deviation = 2*(middle-rejectionratesincelastadjustment)/span;       % positive, if too small rejectionrate
                        if( deviation*obj.adjustedtoofar(j_up)<=0 )                         % i.e. on different side of optimum this time
                            obj.adjustedtoofar(j_up) = abs(obj.adjustedtoofar(j_up)) + 1;   % overadjusted one more time
                        end
                        obj.adjustedtoofar(j_up) = sign(deviation)*abs(obj.adjustedtoofar(j_up));	% switch sign, so it fits with deviation
                        obj.stepsizes(j_up) = obj.stepsizes(j_up)*(adjfactor^(deviation/(2^abs(obj.adjustedtoofar(j_up)))));
                        if( deviation>0 )       % i.e. too low rejectionrate
                            fprintf( ' (%s) Info - optimisestepsizes (%d): Increased step%d to %1.1e (rej. rate = %1.3e)\n', obj.comment,obj.MCit,j_up, obj.stepsizes(j_up), rejectionratesincelastadjustment );
                        elseif( deviation<0 )   % i.e. too high rejectionrate
                            fprintf( ' (%s) Info - optimisestepsizes (%d): Decreased step%d to %1.1e (rej. rate = %1.3e)\n', obj.comment,obj.MCit,j_up, obj.stepsizes(j_up), rejectionratesincelastadjustment );
                        end     % end of checking of rejectionrates reasonable
                        % reset counter:
                        obj.rejected(j_up) = 0; obj.samplecounter(j_up) = 0;
                    end     % end if enough samples
                end     % end of going through updates
                
                % reset timer and memory:
                obj.nextoptstepstime = obj.MCit + optstepsfreq;
            end     % end of checking burnin
            %diary 'getgeometry controlwindow output.txt';
     	end     % end of optimisestepsizes function
        function[] = standardisetruegeometry( obj )
            % standardises the true geometry, so it has fl1 at origin, fl2 at x-axis and fl3 in x-y-plane
            % only updates the true geometry and perspective current values, not the proposed values
            
            % set auxiliary parameters:
            datarange = 1:obj.nodata;           % all data possible
            
            for j_state = 1:obj.nostates
                % standardise true geometry:
                % ...translations to origin:
                overalltrans = -obj.flpos_curr(:,1,j_state);
                obj.flpos_curr(:,:,j_state) = obj.flpos_curr(:,:,j_state) + overalltrans;
                % ...rotation around fl1 so fl2-fl1 is on x-axis:
                fl1fl2 = (obj.flpos_curr(:,2,j_state)-obj.flpos_curr(:,1,j_state));   fl1fl2 = fl1fl2/norm(fl1fl2);
                a = cross( fl1fl2, [1;0;0] ); alpha = acos(fl1fl2'*[1;0;0]);    % no minus sign here
                R3d = obj.get3drotation( a,alpha );
                obj.flpos_curr(:,:,j_state) = R3d*(obj.flpos_curr(:,:,j_state)-obj.flpos_curr(:,1,j_state)) + obj.flpos_curr(:,1,j_state);
                overallR3d = R3d;
                % ...rotation around fl1(s) so fl3-fl1 is in xy-plane:
                if( obj.nomarkers>2 )
                    a = [1;0;0];        % = fl1-fl2-axis
                    projection = (obj.flpos_curr(:,3,j_state)-obj.flpos_curr(:,1,j_state));
                    projection = projection - (projection'*a)*a; projection = projection/norm(projection);
                    alpha = sign((a'*cross(projection,[0;1;0])))*acos(projection'*[0;1;0]); % get rotation sign, from orientation of cross product 
                    R3d = obj.get3drotation( a,alpha );
                    obj.flpos_curr(:,:,j_state) = R3d*(obj.flpos_curr(:,:,j_state)-obj.flpos_curr(:,1,j_state)) + obj.flpos_curr(:,1,j_state);
                    overallR3d = R3d*overallR3d;
                end     % end if more than two fluorophores
                
                % adjust perspectives:
                obj.trans_curr(:,obj.statevskin_curr(j_state,:)) = obj.trans_curr(:,obj.statevskin_curr(j_state,:)) - repmat(overalltrans,[1,sum(obj.statevskin_curr(j_state,:))]);
                for j_data = datarange(obj.statevskin_curr(j_state,:))
                    obj.rot_curr(:,:,j_data) = obj.rot_curr(:,:,j_data)*(overallR3d');   % overallR3d' is inverse of overallR3d
                end     % end of data loop
            end     % end of going through states
    	end     % end of standardisetruegeometry function
        function[] = keepmeansinvariant( obj, j_data )
            % adjusts translations so proposed positions don't change the
            % mean of the predictions, compared to the current positioins
            
            % set auxiliary parameters:
            if( j_data<0 )
                datarange = 1:obj.nodata;
            else
                datarange = j_data;
            end
            
            % get current means:
            mean_curr = mean( obj.flpos_pred_curr(:,:,datarange),2 );   % mean over markers
            % get proposed means:
            obj.updatepositionpredictions(j_data,-1);                   % update proposed predictions
            mean_prop = mean( obj.flpos_pred_prop(:,:,datarange),2 );   % mean over markers
            obj.trans_prop(:,datarange) = obj.trans_curr(:,datarange) - squeeze(mean_prop-mean_curr);
        end     % end of keepmeansinvariant function
        
        function[] = regularcontrolwindowoutput( obj )
            % updates the control window information with current state
            if( mod(obj.MCit,obj.outputfreq)==0 || obj.MCit==obj.burnin || obj.MCit==obj.MCmax )
                fprintf( ' (%s)...current output at MCit=%d (%1.3f sec)(logtarget_curr = %+1.8e):\n', obj.comment,obj.MCit, (now-obj.t1)*24*3600, obj.logtarget_curr );
                for j_state = 1:obj.nostates
                    fprintf( ' (%s) ...state %d:\t(%d measurements, [ %s])\n', obj.comment,j_state, sum(obj.statevskin_curr(j_state,:)), sprintf('%5.2f%% ',100*obj.statefrac_curr(j_state,:)) );
                    % ...lengths:
                    for j_marker = 1:(obj.nomarkers-1)
                        fprintf( ' (%s)  ', obj.comment );
                        for jj_marker = (j_marker+1):obj.nomarkers
                            fprintf( 'l%d%d = %1.3e\t', j_marker,jj_marker, sqrt(sum((obj.flpos_curr(:,j_marker,j_state)-obj.flpos_curr(:,jj_marker,j_state)).^2,1)));
                        end     % end of inner markers loop
                        fprintf('\n');
                    end     % end of outer markers loop
                    fprintf( ' (%s)  ', obj.comment );
                    % ...measurement errors:
                    for j_marker = 1:obj.nomarkers
                        fprintf( 'merr%d = [ %1.2e %1.2e ]\t', j_marker, obj.merr_curr(1,1,j_marker,j_state), obj.merr_curr(3,3,j_marker,j_state) );
                    end     % end of markers loop
                    fprintf('\n');
                end     % end of states loop
                if( obj.nostates>1 )
                    % ...state fractions:
                    for j_dset = 1:obj.nodsets
                        fprintf( ' (%s) statefrac(%d) = [ %s]\n',obj.comment,j_dset, sprintf('%1.2e ',obj.statefrac_curr(:,j_dset)) );
                    end     % end of datasets loop
                    % ...occupation numbers:
                    stateoptions = 1:obj.nostates;
                    word = sprintf('%d ',stateoptions*obj.statevskin_curr);  if( word>100 ), word = [word(1:97),'...']; end
                    fprintf( ' (%s) statevskin  = [ %s]\n', obj.comment, word );
                end     % end if label model multistate
                fprintf( ' (%s) rejection rates: [ %s]\n', obj.comment, sprintf( '%6.4f ', obj.rejected./obj.samplecounter )  );
                %diary 'getgeometry controlwindow output.txt';
            end     % end if output now
        end     % end of regularcontrolwindowoutput function
        function[] = drawthisprediction( obj )
            % draws a ball centred at each mean fluorophore position and
            % with radius equal to the measurement error
            
            depictdata = min(1,obj.nodata);                 % number of datapoints that actually gets depicted
            stateoptions = 1:obj.nostates;                  % all states possible
            
            fig = figure('Name','Current predictions');
            title( sprintf('current sample for predictions (%s)\n(MCit=%d)(merr=[%1.3e,%1.3e,%1.3e])', obj.comment, obj.MCit, sprintf( '%1.3e ',obj.merr_curr(1,1,:,:) ) ) );
            xlabel('x'); ylabel('y'); zlabel('z'); %axis equal;
            %xlim([-4000,4000]); ylim([-4000,4000]); zlim([-4000,4000]);
            grid on; box on; view([30,30]);
            hold on;
            
            % predictions:
            [x,y,z] = sphere(6);
            for j_data = 1:depictdata
                j_state = stateoptions(obj.statevskin_curr(:,j_data));
                if( ~isempty(j_state) )
                    for j_marker = 1:obj.nomarkers
                        x_here = x*obj.merr_curr(1,1,j_marker,j_state); y_here = y*obj.merr_curr(2,2,j_marker,j_state); z_here = z*obj.merr_curr(3,3,j_marker,j_state);
                        thisx = x_here + obj.flpos_pred_curr(1,j_marker,j_data);    thisy = y_here + obj.flpos_pred_curr(2,j_marker,j_data);    thisz = z_here + obj.flpos_pred_curr(3,j_marker,j_data);
                        surf( thisx,thisy,thisz, 'FaceAlpha',0.1, 'FaceColor','blue', 'EdgeColor','none' );
                    end     % end of markers loop
                end     % end if state exists for this datapoint
            end     % end of going through all data indices
            
            % observations:
            buffer = reshape(obj.flpos_data(:,:,1:depictdata),[3,obj.nomarkers*depictdata]);
            colorscheme = 12000*ceil((1:size(buffer,2))/obj.nodata)-6000;
            scatter3( buffer(1,:), buffer(2,:), buffer(3,:), 20, colorscheme, 'filled', 'o' ); xlabel('x'); ylabel('y'); zlabel('z'); title( 'all fluorophores' );  colormap jet; colorbar;
            for j_data = 1:depictdata
                line( obj.flpos_data(1,[1:obj.nomarkers,1],j_data),obj.flpos_data(2,[1:obj.nomarkers,1],j_data),obj.flpos_data(3,[1:obj.nomarkers,1],j_data), 'Color','black' );
            end     % end of going through kinetochores
            hold off; drawnow;  obj.getprinted( fig, get(get(gca,'title'),'string') );
            
            % for deviations:
            for j_marker = 1:obj.nomarkers
                fig = figure; values = squeeze(obj.flpos_data(:,j_marker,:)-obj.flpos_pred_curr(:,j_marker,:)); scatter3( values(1,:),values(2,:),values(3,:), 15,'filled','o' ); xlabel('x'); ylabel('y'); zlabel('z'); title( sprintf('fl%d deviations (%s)\nmean=[%+1.2f,%+1.2f,%+1.2f],eigstd=[%1.2f,%1.2f,%1.2f]', j_marker,obj.comment, mean(values,2),sqrt(eig(cov(values')))) ); axis equal;
                drawnow; obj.getprinted( fig, get(get(gca,'title'),'string') );
            end     % end of markers loop
        end     % end of drawthisprediction function
        function[] = drawaverageprediction( obj )
            % same as drawthisprediction but with overall average up to this point
            
            % remember current predictions:
            buffer = obj.flpos_pred_curr;   % to be able to memorise
            
            % set current predictions to average predictions, in order to recycle the drawthisprediction function:
            if( obj.MCit<obj.burnin+1 ) % still during burnin
                statsrangehere = 1:(obj.MCit + 0);         % don't include latest one, just in case it's before the update
            else                        % only considre post-burnin data
                statsrangehere = obj.statsrange(1):(obj.MCit + 0);
            end
            for j_it = 1:numel(statsrangehere)
                jj_it = statsrangehere(j_it);
                obj.flpos_prop(:,:,:) = obj.flpos_hist(:,:,:,jj_it);
                obj.trans_prop(:,:) = obj.trans_hist(:,:,jj_it);
                obj.rot_prop(:,:,:) = obj.rot_hist(:,:,:,jj_it);
                obj.updatepositionpredictions(-1,-1);
                % compute floating averages:
                obj.flpos_pred_curr(:,:,:) = ( obj.flpos_pred_curr*(j_it-1) + obj.flpos_pred_prop*1 )/j_it;
            end     % end MCit loop
            
            % draw the 'current' prediction:
            obj.drawthisprediction();
            
            % reset current predictions to original situation:
            obj.flpos_pred_curr = buffer;
            obj.copycurrtoprops();	% reset proposal
        end     % end of drawaverageprediction function
        function[] = outputlengthstatistics( obj, withcw, withhist, withscat, withevol )
            % plots the histograms for the (naive) lengths between the fluorophore positions
            
            for j_state = 1:obj.nostates
                % get lengths:
                l_hist = zeros(obj.nomarkers,obj.nomarkers,obj.MCmax);
                for j_marker = 1:(obj.nomarkers-1)
                    for jj_marker = (j_marker+1):obj.nomarkers
                        l_hist(j_marker,jj_marker,:) = squeeze( sqrt(sum((obj.flpos_hist(:,j_marker,j_state,:)-obj.flpos_hist(:,jj_marker,j_state,:)).^2,1)) );
                        l_hist(jj_marker,j_marker,:) = l_hist(j_marker,jj_marker,:);    % symmetric
                    end     % end of inner markers loop
                end     % end of outer markers loop
                % output:
                if( withcw==1 )         % output in control-window
                    fprintf( '  Info - outputlengthstatistics (%s): Some length-statistics for state %d (after burnin of %d):\n', obj.comment, j_state, obj.burnin );
                    for j_marker = 1:(obj.nomarkers-1)
                        for jj_marker = (j_marker+1):obj.nomarkers
                            fprintf( '   l%d%d  = %1.3e +- %1.3e\t(%s)\n', j_marker,jj_marker, mean(l_hist(j_marker,jj_marker,obj.statsrange)), std(l_hist(j_marker,jj_marker,obj.statsrange)), obj.comment );
                        end     % end of inner markers loop
                    end     % end of outer markers loop
                    %diary 'getgeometry controlwindow output.txt';
                end     % end of withcw
                if( withhist==1 )       % output marginal histograms
                    % ...lengths:
                    res = ceil( 4*numel(obj.statsrange)^(1/3) );
                    for j_marker = 1:(obj.nomarkers-1)
                        for jj_marker = (j_marker+1):obj.nomarkers
                            fig = figure;   title(sprintf('l%d%d samples (%d)(after burnin)',j_marker,jj_marker,j_state));   xlabel(sprintf('l%d%d',j_marker,jj_marker));   ylabel('freq');    hold on;   histogram( l_hist(j_marker,jj_marker,obj.statsrange),res );	hold off;
                            drawnow; obj.getprinted( fig, get(get(gca,'title'),'string') );
                        end     % end of inner markers loop
                    end     % end of outer markers loop
                    % ...measurement errors:
                    for j_marker = 1:obj.nomarkers
                        fig = figure; histogram( obj.merr_hist(1,1,j_marker,j_state,obj.statsrange), res );     title(sprintf('merr%dxy samples (%d)(after burnin)',j_marker,j_state));	xlabel(sprintf('merr%dxy',j_marker));   ylabel('freq');
                        drawnow; obj.getprinted( fig, get(get(gca,'title'),'string') );
                        fig = figure; histogram( obj.merr_hist(3,3,j_marker,j_state,obj.statsrange), res );     title(sprintf('merr%dz  samples (%d)(after burnin)',j_marker,j_state));	xlabel(sprintf('merr%dz',j_marker));    ylabel('freq');
                        drawnow; obj.getprinted( fig, get(get(gca,'title'),'string') );
                    end     % end of markers loop
                    if( obj.nostates>1 )
                        for j_dset = 1:obj.nodsets
                            fig = figure; histogram( obj.statefrac_hist(j_state,j_dset,obj.statsrange), res );     title( sprintf('state %d fraction in set %d (after burnin)',j_state,j_dset) );     xlabel(sprintf('state %d frac in set %d',j_state,j_dset));     ylabel('freq');
                            drawnow; obj.getprinted( fig, get(get(gca,'title'),'string') );
                        end     % end of datasets loop
                    end     % end if multistates model
                end     % end of withhist
                if( (withscat==1) && (obj.nomarkers>2) )	% output scatter plot of the three triangle lengths
                    coval = sqrt(sum(squeeze( obj.merr_hist(1,1,:,obj.statsrange).^2 + obj.merr_hist(2,2,:,obj.statsrange).^2 + obj.merr_hist(3,3,:,obj.statsrange).^2 ),1));  % sum over markers
                    fig = figure; scatter3( l_hist(1,2,obj.statsrange),l_hist(1,3,obj.statsrange),l_hist(2,3,obj.statsrange), 5,coval, 'filled','o' ); title(sprintf('triangle lengths (%d)(after burnin)',j_state)); xlabel('l12'); ylabel('l13'); zlabel('l23'); h = colorbar; set(get(h,'title'),'string','errsum');
                    drawnow; obj.getprinted( fig, get(get(gca,'title'),'string') );
                end     % end of withscat
                if( withevol==1 )       % output of evolution plots of MC
                    MCindices = 1:obj.MCmax;
                    % ...lengths:
                    fig = figure; title(sprintf('MC evolution of lengths (%d)',j_state)); xlabel('MC iteration'); ylabel('length'); hold on;
                    for j_marker = 1:(obj.nomarkers-1)
                        for jj_marker = (j_marker+1):obj.nomarkers
                            plot(MCindices,squeeze(l_hist(j_marker,jj_marker,MCindices)));
                            sofar = get(gca,'legend');
                            if( isempty(sofar)==0 )
                                legend( sofar.String{1:(end-1)}, sprintf('l%d%d',j_marker,jj_marker) );
                            else
                                legend( sprintf('l%d%d',j_marker,jj_marker) );
                            end
                        end     % end of inner markers loop
                    end     % end of outer markers loop
                    hold off; drawnow; obj.getprinted( fig, get(get(gca,'title'),'string') );
                    % ...measurement errors:
                    for j_marker = 1:obj.nomarkers
                        fig = figure; title(sprintf('MC evolution of measurement errors (%d)',j_state)); xlabel('MC iteration'); ylabel(sprintf('merr%d',j_marker)); hold on;
                        plot(MCindices,squeeze(obj.merr_hist(1,1,j_marker,j_state,MCindices)));
                        plot(MCindices,squeeze(obj.merr_hist(3,3,j_marker,j_state,MCindices)));
                        legend(sprintf('merr%dxy',j_marker),sprintf('merr%dz',j_marker));
                        hold off;  drawnow; obj.getprinted( fig, get(get(gca,'title'),'string') );
                    end     % end of markers loop
                    % ...state fractions:
                    if( (j_state==1) && (obj.nostates>1) )  % first of multiple states
                        for j_dset = 1:obj.nodsets
                            fig = figure; title( sprintf('state fractions, set %d',j_dset) ); xlabel('MC iteration'); ylabel('state frac'); hold on;
                            for jj_state = 1:obj.nostates
                                plot( MCindices,obj.statefrac_hist(jj_state,j_dset,MCindices) );
                                sofar = get(gca,'legend');
                                if( isempty(sofar)==0 )
                                    legend( sofar.String{1:(end-1)}, sprintf('state %d',jj_state) );
                                else
                                    legend( sprintf('state %d',jj_state) );
                                end
                            end     % end of label states loop
                            hold off; drawnow; obj.getprinted( fig, get(get(gca,'title'),'string') );
                        end     % end of datasets loop
                    end     % end if multistate model
                end     % end of withevol
            end     % end of state loop
        end     % end of outputlengthstatistics function
        function[] = outputinternalanglesstatistics( obj, withcw, withhist, withscat, withevol )
            % plots histograms of the inferred, true and naive internal angles
            
            % set auxiliary parameters:
            if( obj.nomarkers~=3 ), return; end     % only meant to work for triangles
            outrange = 1:obj.MCmax;	% outputrange for evolution
            res = ceil(0.5*(obj.nodata*numel(obj.statsrange))^(1/3));
            
            for j_state = 1:obj.nostates
                % reproduce internal angles:
                % ...angle at fl1:
                r1 = squeeze(obj.flpos_hist(:,2,j_state,outrange)-obj.flpos_hist(:,1,j_state,outrange)); r1 = r1./sqrt(sum(r1.*r1));
                r2 = squeeze(obj.flpos_hist(:,3,j_state,outrange)-obj.flpos_hist(:,2,j_state,outrange)); r2 = r2./sqrt(sum(r2.*r2));
                phi213_pred = acos(sum(r1.*r2));                                % inferred angles at fl1
                r1 = squeeze(obj.flpos_data(:,2,:)-obj.flpos_data(:,1,:)); r1 = r1./sqrt(sum(r1.*r1));
                r2 = squeeze(obj.flpos_data(:,3,:)-obj.flpos_data(:,1,:)); r2 = r2./sqrt(sum(r2.*r2));
                phi213_naive_prel = acos(sum(r1.*r2));                      	% naive angles at fl1
                % ...angle at fl2:
                r1 = squeeze(obj.flpos_hist(:,3,j_state,outrange)-obj.flpos_hist(:,2,j_state,outrange)); r1 = r1./sqrt(sum(r1.*r1));
                r2 = squeeze(obj.flpos_hist(:,1,j_state,outrange)-obj.flpos_hist(:,2,j_state,outrange)); r2 = r2./sqrt(sum(r2.*r2));
                phi321_pred = acos(sum(r1.*r2));                                % inferred angles at fl2
                r1 = squeeze(obj.flpos_data(:,3,:)-obj.flpos_data(:,2,:)); r1 = r1./sqrt(sum(r1.*r1));
                r2 = squeeze(obj.flpos_data(:,1,:)-obj.flpos_data(:,2,:)); r2 = r2./sqrt(sum(r2.*r2));
                phi321_naive_prel = acos(sum(r1.*r2));                         	% naive angles at fl2
                % ...angle at fl3:
                r1 = squeeze(obj.flpos_hist(:,1,j_state,outrange)-obj.flpos_hist(:,3,j_state,outrange)); r1 = r1./sqrt(sum(r1.*r1));
                r2 = squeeze(obj.flpos_hist(:,2,j_state,outrange)-obj.flpos_hist(:,3,j_state,outrange)); r2 = r2./sqrt(sum(r2.*r2));
                phi132_pred = acos(sum(r1.*r2));                                % inferred angles at fl3
                r1 = squeeze(obj.flpos_data(:,1,:)-obj.flpos_data(:,3,:)); r1 = r1./sqrt(sum(r1.*r1));
                r2 = squeeze(obj.flpos_data(:,2,:)-obj.flpos_data(:,3,:)); r2 = r2./sqrt(sum(r2.*r2));
                phi132_naive_prel = acos(sum(r1.*r2));                         	% naive angles at fl3

                % controlwindow output:
                if( withcw==1 )
                    fprintf( '  Info - outputinternalanglesstatistics (%s): Some internal angles-statistis for state %d (after burnin of %d)\n', obj.comment, j_state, obj.burnin );
                    fprintf( '   phi213(%d) = %1.3e +- %1.3e\t(%s)\n', j_state, mean(phi213_pred(obj.statsrange))*180/pi, std(phi213_pred(obj.statsrange))*180/pi, obj.comment );
                    fprintf( '   phi321(%d) = %1.3e +- %1.3e\t(%s)\n', j_state, mean(phi321_pred(obj.statsrange))*180/pi, std(phi321_pred(obj.statsrange))*180/pi, obj.comment );
                    fprintf( '   phi132(%d) = %1.3e +- %1.3e\t(%s)\n', j_state, mean(phi132_pred(obj.statsrange))*180/pi, std(phi132_pred(obj.statsrange))*180/pi, obj.comment );
                    %diary 'getgeometry controlwindow output.txt';
                end     % end if withcw
                % graphical output:
                % ...of histograms:
                if( withhist==1 )
                    % ....sample from naive preliminary (i.e. true) values, so equal amount as for inferred samples:
                    phi213_naive = zeros(1,length(obj.statsrange)); phi321_naive = zeros(1,length(obj.statsrange)); phi132_naive = zeros(1,length(obj.statsrange));
                    for j_sample = 1:numel(obj.statsrange)
                        index = ceil(rand()*obj.nodata);
                        phi213_naive(j_sample) = phi213_naive_prel(index);
                        phi321_naive(j_sample) = phi321_naive_prel(index);
                        phi132_naive(j_sample) = phi132_naive_prel(index);
                    end
                    % ....plot:
                    fig = figure;   title( sprintf('phi213 (%d)(after burnin)',j_state) );  xlabel('angle in deg'); hold on; histogram( phi213_naive*180/pi,res, 'FaceColor','yellow');   histogram( phi213_pred(obj.statsrange)*180/pi,res,'FaceColor','blue');    hold off;  legend('naive','pred');
                    drawnow; obj.getprinted( fig, get(get(gca,'title'),'string') );
                    fig = figure;   title( sprintf('phi321 (%d)(after burnin)', j_state) );	xlabel('angle in deg'); hold on; histogram( phi321_naive*180/pi,res, 'FaceColor','yellow');   histogram( phi321_pred(obj.statsrange)*180/pi,res,'FaceColor','blue');    hold off;  legend('naive','pred');
                    drawnow; obj.getprinted( fig, get(get(gca,'title'),'string') );
                    fig = figure;   title( sprintf('phi132 (%d)(after burnin)', j_state) ); xlabel('angle in deg'); hold on; histogram( phi132_naive*180/pi,res, 'FaceColor','yellow');   histogram( phi132_pred(obj.statsrange)*180/pi,res,'FaceColor','blue');    hold off;  legend('naive','pred');
                    drawnow; obj.getprinted( fig, get(get(gca,'title'),'string') );
                end     % end if withhist
                % ...of scatterplots:
                if( withscat==1 )
                    fig = figure; scatter3( phi213_pred(obj.statsrange)*180/pi, phi321_pred(obj.statsrange)*180/pi, phi132_pred(obj.statsrange)*180/pi, 5, 'filled', 'o' ); title( sprintf('phi213 vs phi321 vs phi132 (%d)(after burnin)',j_state) ); xlabel('phi213'); ylabel('phi321'); zlabel('phi132');
                    drawnow; obj.getprinted( fig, get(get(gca,'title'),'string') );
                end     % end if withscat
                % ...of evolutionplots:
                if( withevol==1 )
                    MCindices = 1:obj.MCmax;
                    fig = figure; plot( MCindices, phi213_pred(outrange)*180/pi ); hold on; plot( MCindices, phi321_pred(outrange)*180/pi ); plot( MCindices, phi132_pred(outrange)*180/pi ); hold off; title( sprintf('internal angles evolution (%d)',j_state) ); legend('phi213','phi321','phi132'); xlabel('MCit');
                    drawnow; obj.getprinted( fig, get(get(gca,'title'),'string') );
                end     % end if withevol
            end     % end of states loop
        end     % end of outputinternalanglesstatistics function
        function[] = outputcurrentandproposedpars( obj )
            % outputs current states and proposals for true geometry,
            % error, perspectives and predictions
            
            fprintf( ' Info - outputcurrentandproposedpars (%s)(%d):\n', obj.comment,obj.MCit );
            % logtarget:
            fprintf( '   logtarget_prop = %+1.5e\n', obj.logtarget_prop );
            fprintf( '   logtarget_curr = %+1.5e\n', obj.logtarget_curr );
            for j_state = 1:obj.nostates
                % true geometry:
                fprintf( '  template positions:\n' );
                for j_marker = 1:obj.nomarkers
                    fprintf('   fl%d_curr(%d) = [ %s]\n', j_marker,j_state, sprintf('%+1.3e',obj.flpos_curr(:,j_marker,j_state)) );
                    fprintf('   fl%d_prop(%d) = [ %s]\n', j_marker,j_state, sprintf('%+1.3e',obj.flpos_prop(:,j_marker,j_state)) );
                end     % end of markers loop
                for j_marker = 1:(obj.nomarkers-1)
                    for jj_marker = (j_marker+1):obj.nomarkers
                        fprintf( '   l%d%d_curr(%d) = %1.3e\n', j_marker,jj_marker,j_state, sqrt(sum((obj.flpos_curr(:,j_marker,j_state)-obj.flpos_curr(:,jj_marker,j_state)).^2,1)) );
                        fprintf( '   l%d%d_prop(%d) = %1.3e\n', j_marker,jj_marker,j_state, sqrt(sum((obj.flpos_prop(:,j_marker,j_state)-obj.flpos_prop(:,jj_marker,j_state)).^2,1)) );
                    end     % end of inner markers loop
                end     % end of outer markers loop
                % error:
                fprintf( '  error:\n' );
                for j_marker = 1:obj.nomarkers
                    fprintf( '   merr%d_curr(%d) = [ %1.3e %1.3e ]\n', j_marker,j_state, obj.merr_curr(1,1,j_marker,j_state),obj.merr_curr(3,3,j_marker,j_state) );
                    fprintf( '   merr%d_prop(%d) = [ %1.3e %1.3e ]\n', j_marker,j_state, obj.merr_prop(1,1,j_marker,j_state),obj.merr_prop(3,3,j_marker,j_state) );
                end     % end of markers loop
            end     % end of states loop
            % perspectives:
            datarange = [1,ceil(obj.nodata/2)];
            fprintf( '  perspectives and predictions:\n' );
            for j_data = datarange
                fprintf( '   trans_curr%d = [ %s]\n', j_data, sprintf('%+1.3e ',obj.trans_curr(:,j_data)) );
                fprintf( '   trans_prop%d = [ %s]\n', j_data, sprintf('%+1.3e ',obj.trans_prop(:,j_data)) );
                fprintf( '   rot_curr%d   = [ %s]\n', j_data, sprintf('%+1.3e ',obj.rot_curr(:,:,j_data)) );
                fprintf( '   rot_prop%d   = [ %s]\n', j_data, sprintf('%+1.3e ',obj.rot_prop(:,:,j_data)) );
                for j_marker = 1:obj.nomarkers
                    fprintf('   fl%d_pred_curr%d = [ %s]\n', j_marker,j_data, sprintf('%+1.3e ',obj.flpos_pred_curr(:,j_marker,j_data)) );
                    fprintf('   fl%d_pred_prop%d = [ %s]\n', j_marker,j_data, sprintf('%+1.3e ',obj.flpos_pred_prop(:,j_marker,j_data)) );
                end     % end of markers loop
            end     % end of data loop
            %diary 'getgeometry controlwindow output.txt';
            pause(2);
        end     % end of outputcurrentandproposedpars function
        function[] = comparecurrentandproposedpars( obj )
            % outputs current states and proposals for true geometry,
            % error, perspectives and predictions
            
            cond = (0==1);      % turns '1' as soon as one difference detected
            % logtarget:
            %diff = max(abs(obj.logtarget_prop-obj.logtarget_curr),[],'all');
            %if( diff>0 ), fprintf(' Info - comparecurrentandproposedpars (%s)(%d): logtarget different, diff = %+1.7e\n', obj.comment,obj.MCit, diff); cond = (1==1); end
            for j_state = 1:obj.nostates
                % true geometry:
                for j_marker = 1:obj.nomarkers
                    diff = max(abs(obj.flpos_prop(:,j_marker,j_state)-obj.flpos_curr(:,j_marker,j_state)),[],'all');
                    if( diff>0 ), fprintf(' Info - comparecurrentandproposedpars (%s)(%d): template %d,%d different, diff = %+1.7e\n', obj.comment,obj.MCit, j_marker,j_state, diff); cond = (1==1); end
                end     % end of markers loop
                % error:
                for j_marker = 1:obj.nomarkers
                    diff = max(abs(obj.merr_prop(:,:,j_marker,j_state)-obj.merr_curr(:,:,j_marker,j_state)),[],'all');
                    if( diff>0 ), fprintf(' Info - comparecurrentandproposedpars (%s)(%d): error %d,%d different, diff = %+1.7e\n', obj.comment,obj.MCit, j_marker,j_state, diff); cond = (1==1); end
                end     % end of markers loop
            end     % end of states loop
            % perspectives:
            datarange = 1:obj.nodata;
            for j_data = datarange
                diff = max(abs(obj.trans_prop(:,j_data)-obj.trans_curr(:,j_data)),[],'all');
                if( diff>0 ), fprintf(' Info - comparecurrentandproposedpars (%s)(%d): trans %d different, diff = %+1.7e\n', obj.comment,obj.MCit, j_data, diff); cond = (1==1); end
                diff = max(abs(obj.rot_prop(:,:,j_data)-obj.rot_curr(:,:,j_data)),[],'all');
                if( diff>0 ), fprintf(' Info - comparecurrentandproposedpars (%s)(%d): rot %d different, diff = %+1.7e\n', obj.comment,obj.MCit, j_data, diff); cond = (1==1); end
                for j_marker = 1:obj.nomarkers
                    diff = max(abs(obj.flpos_pred_prop(:,j_marker,j_data)-obj.flpos_pred_curr(:,j_marker,j_data)),[],'all');
                    if( diff>0 ), fprintf(' Info - comparecurrentandproposedpars (%s)(%d): pred %d,%d different, diff = %+1.7e\n', obj.comment,obj.MCit, j_marker,j_data, diff); cond = (1==1); end
                end     % end of markers loop
            end     % end of data loop
            % states:
            if( obj.nostates>1 )
                diff = max(abs(obj.statefrac_prop-obj.statefrac_curr),[],'all');
                if( diff>0 ), fprintf(' Info - comparecurrentandproposedpars (%s)(%d): statefrac different, diff = %+1.7e\n', obj.comment,obj.MCit, diff); cond = (1==1); end
                diff = max(abs(obj.statevskin_prop-obj.statevskin_curr),[],'all');
                if( diff>0 ), fprintf(' Info - comparecurrentandproposedpars (%s)(%d): statevskin different, diff = %+1.7e\n', obj.comment,obj.MCit, diff); cond = (1==1); end
            end     % end if multiple states
            if( cond ), obj.outputcurrentandproposedpars(); end
            %diary 'getgeometry controlwindow output.txt';
        end     % end of outputcurrentandproposedpars function
        function[] = outputmarginalstates( obj, withcw,withhist )
            % outputs the marginal distribution over the states for each datapoint
            
            if( obj.nostates>1 )   % i.e. multistate model
                state_histogram = squeeze(mean(obj.statevskin_hist(:,:,obj.statsrange),3));
                sttprops = squeeze( mean( obj.statevskin_hist(:,:,obj.statsrange),2 ) );
                
                if( withcw==1 )
                    fprintf( ' Info - outputmarginalstates (%d)(%s): state marginals per datapoint (after burnin of %d) are:\n', obj.MCit,obj.comment, obj.burnin );
                    for j_state = 1:obj.nostates
                        fprintf( '  %d:\t%1.3e +- %1.3e,\t[ %s]\t(%s)\n', j_state, mean(sttprops(j_state,:),2), std(sttprops(j_state,:),0,2), sprintf('%1.2f ',state_histogram(j_state,:)), obj.comment );
                    end     % end of states loop
                end     % end of withcw
                if( withhist==1 )
                    fig = figure; title( sprintf( 'state marginals per datapoint (after burnin)' ) ); colorbar; yticks(1:obj.nostates); hold on;
                    image( state_histogram, 'CDataMapping','scaled' ); 
                    hold off; drawnow; obj.getprinted( fig, get(get(gca,'title'),'string') );
                end     % end of withhist
            end     % end of distinguishing mutlistate model types
        end     % end of outputmarginalstates function
        function[] = writestatestext( obj )
            % writes current state into text file
            
            textfile = fopen( sprintf('%s/%d-%02.0f-%02.0f_%02.0f-%02.0f_%s.txt',obj.subfoldername, obj.globaltimestamp(1:5),obj.comment), 'a' );
            if( obj.MCit==1 )   % also write header
                fprintf(textfile, '  version:\t\t\t2\n' );
                fprintf(textfile, '%4d-%02d-%02d--%02d-%02d\n', obj.globaltimestamp(1:5) );
                fprintf(textfile, '  filename:\t\t\t%s\n', obj.inputfilename );
                fprintf(textfile, '  comment:\t\t\t%s\n', '' );
                fprintf(textfile, '  chaincomment:\t\t\t%s\n', obj.comment );
                fprintf(textfile, '  nostates:\t\t\t%d\n', obj.nostates );
                fprintf(textfile, '  nomarkers:\t\t\t%d\n', obj.nomarkers );
                fprintf(textfile, '  dim:\t\t\t\t%d\n', 3 );
                fprintf(textfile, '  nodatasets:\t\t\t%d\n', obj.nodsets );
                fprintf(textfile, '  nodata:\t\t\t%d\n', obj.nodata );
                fprintf(textfile, '  samplerange:\t\t[%d..%d]\n', obj.burnin,obj.MCmax );
                fprintf(textfile, '  posstatevskin:\n' );
                for j_state = 1:obj.nostates
                    fprintf(textfile, '%d\t', obj.posstatevskin(j_state,:) );
                    fprintf(textfile, '\n' );
                end     % end of nostates
                fprintf(textfile, '  datasetvskin:\n' );
                for j_dset = 1:obj.nodsets
                    fprintf(textfile, '%d\t', obj.datasetvskin(j_dset,:) );
                    fprintf(textfile, '\n' );
                end     % end of nostates
                fprintf(textfile, '  lengthprior:\trectangle\n' );
                for j_state = 1:obj.nostates
                    fprintf(textfile, '%+1.5e\t', obj.lengthpriordensity.pars{j_state}' );
                    fprintf(textfile, '\n' );
                end     % end of states loop
                fprintf(textfile, '  merrprior:\tgamma\n' );
                for j_state = 1:obj.nostates
                    for j_marker = 1:obj.nomarkers
                        merrpriors_here = [obj.shapepar_prior{j_state}(j_marker),obj.merrmeans_prior{j_state}(j_marker,:)];
                        fprintf(textfile, '%+1.5e\t', merrpriors_here );
                    end     % end of markers loop
                    fprintf(textfile, '\n' );
                end     % end of states loop
                fprintf(textfile, '  pprior:\t\tnone\n' );
                fprintf(textfile, '%+1.5e\t', zeros(obj.nostates,1) );
                fprintf(textfile, '\n' );
                fprintf(textfile, '\n' );   % five blank lines
                fprintf(textfile, '\n' );
                fprintf(textfile, '\n' );
                fprintf(textfile, '\n' );
                fprintf(textfile, '\n' );
            end     % end if first iteration
            fprintf(textfile, '%d\n', obj.MCit );
            % positions:
            for j_state = 1:obj.nostates
                for j_marker = 1:obj.nomarkers
                    fprintf(textfile, '%+1.5e\t', obj.flpos_curr(:,j_marker,j_state) ); fprintf(textfile, '\n' );
                end     % end of markers loop
            end     % end of states loop
            % measurement errors:
            for j_state = 1:obj.nostates
                for j_marker = 1:obj.nomarkers
                    fprintf(textfile, '%+1.5e\t', diag(obj.merr_curr(:,:,j_marker,j_state)) ); fprintf(textfile, '\n' );
                end     % end of markers loop
            end     % end of states loop
            % translations:
            for j_state = 1:obj.nostates    % other states are duplicates
                for j_data = 1:obj.nodata
                    fprintf(textfile, '%+1.5e\t', obj.trans_curr(:,j_data) ); fprintf(textfile, '\n' );
                end     % end of markers loop
            end     % end of states loop
            % rotations:
            for j_state = 1:obj.nostates    % other states are duplicates
                for j_data = 1:obj.nodata
                    fprintf(textfile, '%+1.5e\t%+1.5e\t%+1.5e\n', obj.rot_curr(:,:,j_data) );
                end     % end of markers loop
            end     % end of states loop
            % p:
            fprintf(textfile, '%+1.5e\t', obj.statefrac_curr' ); fprintf(textfile, '\n' );  % in order st1_dset1, st1_dset2, ..., st2_dset1, ...
            % statevskin:
            fprintf(textfile, '%d\t', (1:obj.nostates)*obj.statevskin_curr ); fprintf(textfile, '\n' );
            % logtarget:
            fprintf(textfile, '%+1.10e\t', obj.logtarget_curr ); fprintf(textfile, '\n' );
            % buffer:
            fprintf(textfile, '\n' );
            fprintf(textfile, '\n' );
            fclose(textfile);
        end     % end of writestatetexts function
    end     % end of private methods
end     % end of getgeometry_v2 class