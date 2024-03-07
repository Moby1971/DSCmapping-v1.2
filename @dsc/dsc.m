classdef dsc


    % Data and parameter class for DSC app
    %
    % Gustav Strijkers
    % g.j.strijkers@amsterdamumc.nl
    % March 2024
    %



    properties

        % Data parameters
        data                        % raw image data
        mask                        % image mask
        dcmInfo                     % Dicom info (tags)
        spatialInfo                 % Dicom spatial info
        position                    % Spatial positions
        positionIndx                % Spatial positions indices
        dataAIF = []                % AIF
        timeAIF = []                % AIF time axes
        ellipseROI = []             % Ellipsoid search region for AIF


        % Dimensions
        dimx = 64                   % nr points readout
        dimy = 64                   % nt points phase encoding
        ns = 1                      % number of slices
        nr = 1                      % number of dynamics
        sl = 1                      % slice thickness
        tr = 500                    % repetition time (ms)
        te = 15                     % echo time (ms)
        fa = 25                     % flip angle
        aspectRatio = 1             % image aspectratio
        pixelAspectRatio = 1        % pixel aspectratio
        fov = 30                    % field-of-view
        maxImageValue = 32767       % max image value
        gifExportSize = 128         % number of horizontal pixels for GIF export


        % Flags
        validDataFlag = false       % valid data file true/false
        validFitFlag = false        % valid DSC fit / calculation


        % Perfusion parameters
        cbv                         % cerebral blood volume
        cbf                         % cerebral blood flow
        mtt                         % mean transit time
        cbvlc                       % cbv with leakage correction
        ttp                         % time to peak
        aif                         % arterial input function
        conc                        % concentration
        s0                          % S0 maps
        fwhm                        % full width half maximum
        k2map                       % leakage parameter


        % Identifiers
        tag                         % file tag for export of dicom files


        % Suppress some messages
        %#ok<*SEPEX>



    end



    methods (Access = public)

        % ---------------------------------------------------------------------------------
        % Empty object constructor
        % ---------------------------------------------------------------------------------
        function obj = dsc()

        end


        % ---------------------------------------------------------------------------------
        % Read the raw dicom data
        % ---------------------------------------------------------------------------------
        function obj = readData(obj, app)

            dcmImportPath = app.dicomImportPath;
            obj.validDataFlag = false;
            obj.validFitFlag = false;

            % List of dicom file names
            flist = dir(fullfile(dcmImportPath,'*.dcm'));

            % Check if dicom data is found
            if ~isempty(flist)

                try

                    app.TextMessage('Reading DICOM images and tags ...');

                    obj.dcmInfo = {};
                    obj.data = [];

                    % Read dicom images
                    [obj.data, obj.spatialInfo, ~] = dicomreadVolume(dcmImportPath);

                    % Read dicom info
                    for indx = 1:length(flist)
                        obj.dcmInfo{indx} = dicominfo([flist(indx).folder,filesep,flist(indx).name]);
                        obj.data(indx) = obj.data(indx)*obj.dcmInfo{indx}.RescaleSlope + obj.dcmInfo{indx}.RescaleIntercept;
                    end

                    % File tag
                    obj.tag = num2str(obj.dcmInfo{1}.SeriesNumber);

                    % Number of unique spatial positions = number of slices
                    [obj.position,~,obj.positionIndx] = unique(obj.spatialInfo.PatientPositions,"rows");
                    nrSlices = size(obj.position,1);

                    % Reshape to dynamics, slices, x, y
                    obj.data = double(obj.data);
                    obj.data = reshape(obj.data,obj.spatialInfo.ImageSize(1),obj.spatialInfo.ImageSize(2),[],nrSlices);
                    obj.data = permute(obj.data,[3,4,2,1]);
                    obj.data = flip(obj.data,4);

                    % Reshape spatial position indices in the same way as the data
                    obj.positionIndx = reshape(obj.positionIndx,1,1,[],nrSlices);
                    obj.positionIndx = permute(obj.positionIndx,[3,4,2,1]);
                    obj.positionIndx = flip(obj.positionIndx,4);
                 
                    % Final image dimensions
                    obj.nr = size(obj.data,1);
                    obj.ns = size(obj.data,2);
                    obj.dimx = size(obj.data,3);
                    obj.dimy = size(obj.data,4);

                    % Aspect ratio
                    obj.pixelAspectRatio = double(obj.dcmInfo{1}.FieldOfViewDimensions(1)*obj.dimy) / double(obj.dcmInfo{1}.FieldOfViewDimensions(2)*obj.dimx);
                    obj.aspectRatio = double(obj.dcmInfo{1}.FieldOfViewDimensions(1)) / double(obj.dcmInfo{1}.FieldOfViewDimensions(2));

                    % Field-of-view
                    obj.fov = obj.dcmInfo{1}.FieldOfViewDimensions(1);

                    % Some sequence parameters
                    obj.sl = obj.dcmInfo{1}.FieldOfViewDimensions(3);
                    obj.tr = obj.dcmInfo{1}.RepetitionTime;
                    obj.te = obj.dcmInfo{1}.EchoTime;
                    obj.fa = obj.dcmInfo{1}.FlipAngle;

                    % Normalize to convenient range
                    obj.data = round(obj.maxImageValue*obj.data/max(obj.data(:)));

                    % AIF
                    obj.dataAIF = zeros(obj.nr,1);
                    obj.timeAIF = (0:obj.nr-1)*obj.tr;  % ms

                    % So far valid data
                    obj.validDataFlag = true;

                    % Some further checks on data validity
                    if size(obj.data,1) < 50
                        app.TextMessage('DSC analysis requires at least 50 dynamics ...');
                        obj.validDataFlag = false;
                    end

                catch ME

                    app.TextMessage(ME.message);

                end

            end

        end % readData



        % ---------------------------------------------------------------------------------
        % Apply threshold to images
        % ---------------------------------------------------------------------------------
        function obj = maskImages(obj, app)

            try

                obj.mask = zeros(obj.ns,obj.dimx,obj.dimy);

                for slice = 1:obj.ns
                    for x = 1:obj.dimx
                        for y = 1:obj.dimy
                            if app.images(1,slice,x,y) > app.threshold(slice)
                                obj.mask(slice,x,y) = 1;
                            end
                        end
                    end
                end

            catch ME

                app.TextMessage(ME.message);

            end

        end % maskImages



        % ---------------------------------------------------------------------------------
        % Update progress bar
        % ---------------------------------------------------------------------------------
        function obj = updateProgressBar(obj, app)

            % Update progress bar
            app.FitProgressGauge.Value = round(100 * app.counter/app.totalNrSteps);
            drawnow;

            % Update the timing indicator
            estimatedtotaltime = app.elapsedTime * app.totalNrSteps / app.counter;
            timeRemaining = estimatedtotaltime * (app.totalNrSteps - app.counter) / app.totalNrSteps;
            timeRemaining(timeRemaining<0) = 0;
            app.EstimatedFitTimeViewField.Value = strcat(datestr(seconds(timeRemaining),'MM:SS')," min:sec"); %#ok<*DATST>
            drawnow;

        end



        % ---------------------------------------------------------------------------------
        % Automatic determination of AIF
        % ---------------------------------------------------------------------------------
        function obj = AIFauto(obj, app)

            % Automatically determines AIF

            % Images: rows, columns, slices, samples
            volumes = permute(app.images,[3,4,2,1]);
            [nR,nC,nS,nT] = size(volumes);

            options = DSC_mri_getOptions(app);
            options.display = 0;
            options.nR = nR;        % Rows
            options.nC = nC;        % Columns
            options.nS = nS;        % Slices
            options.nT = nT;        % Samples (time)

            options.te = obj.te/1000;     % in seconds
            options.tr = obj.tr/1000;     % in seconds
            options.time = (0:nT-1)*obj.tr/1000;  % in seconds
            options.conc = 0;
            options.waitbar = 0;

            % Slice for selecting the AIF
            options.aif.nSlice = app.SlicesEditField.Value;

            % Mask calculation
            autoMask = DSC_mri_mask(volumes, options);

            % Calculate concentrations and S0 maps
            [concMap,~,~] = DSC_mri_conc(volumes, autoMask.data, options);

            % Find AIF
            app.TextMessage("Searching AIF ...");
            options.c1 = app.XViewField.Value;
            options.c2 = app.YViewField.Value;
            [aifResult, roiResult] = DSC_mri_aif(concMap, autoMask.aif, options);

            % Output
            obj.aif = aifResult.conc;
            obj.ellipseROI = roiResult;

            app.roiPoints = [];
            app.roiPoints(:,1) = aifResult.voxels(:,1);
            app.roiPoints(:,2) = aifResult.voxels(:,2);
            app.roiPoints(:,3) = repmat(options.aif.nSlice,length(aifResult.voxels(:,1)),1);

        end



        % ---------------------------------------------------------------------------------
        % Fitting of the DSC indices
        % ---------------------------------------------------------------------------------
        function obj = DSCfit(obj, app)

            % Calculates DSC indices

            % rows, columns, slices, samples
            volumes = permute(app.images,[3,4,2,1]);
            [nR,nC,nS,nT] = size(volumes);

            options = DSC_mri_getOptions(app);
            options.display = 0;
            options.nR = nR;        % Rows
            options.nC = nC;        % Columns
            options.nS = nS;        % Slices
            options.nT = nT;        % Samples (time)

            options.te = obj.te/1000;                     % in seconds
            options.tr = obj.tr/1000;                     % in seconds
            options.time = (0:nT-1)*obj.tr/1000;          % in seconds
            options.conc = 0;
            options.waitbar = 0;

            % Slice for selecting the AIF
            options.aif.nSlice = app.SlicesEditField.Value;

            % Progress bar settings
            app.totalNrSteps = options.step_cvb*options.nS  + options.step_cbv_lc*options.nS + options.step_cbf*options.nS*options.nR + 5;
            app.elapsedTime = 0;
            app.counter = 0;

            % Mask calculation
            app.TextMessage("Masking data ...");
            autoMask = DSC_mri_mask(volumes, options);

            % Calculate concentrations and S0 maps
            app.TextMessage("Calculating concentrations ...");
            [concMap,s0Map,bolus] = DSC_mri_conc(volumes, autoMask.data, options);

            % Find AIF
            options.c1 = app.XViewField.Value;
            options.c2 = app.YViewField.Value;
            [aifResult, roiResult] = DSC_mri_aif(concMap, autoMask.aif, options);

            % CBV calculation
            app.TextMessage("Calculating rCBV ...");
            cbvFit = DSC_mri_cbv(concMap,aifResult.fit.gv, autoMask, options, app);
            obj.cbv = permute(cbvFit,[3,1,2]);

            % CBV leackage correction
            app.TextMessage("Calculating rCBVlc and K2 ...");
            [cbvlcFit,~,K2mapFit,~,~] = DSC_mri_cbv_lc(concMap, aifResult.fit.gv, autoMask, bolus, options, app);
            obj.cbvlc = permute(cbvlcFit,[3,1,2]);
            obj.k2map = permute(K2mapFit,[3,1,2]);

            % CBF calculation
            app.TextMessage("Calculating rCBF ...");
            if app.SVDButton.Value == 1  app.TextMessage('Convolution method = SVD ...');  end
            if app.cSVDButton.Value == 1 app.TextMessage('Convolution method = cSVD ...'); end
            if app.oSVDButton.Value == 1 app.TextMessage('Convolution method = oSVD ...'); end
            cbfFit = DSC_mri_cbf(concMap, aifResult.fit.gv, autoMask, options, app);

            % MTT calculation
            app.TextMessage("Calculating MTT ...");
            mttFit = DSC_mri_mtt(cbvlcFit, cbfFit, options);

            % TTP calculation
            app.TextMessage("Calculating TTP ...");
            ttpFit = DSC_mri_ttp(concMap, autoMask.data, options);
            fwhmFit = DSC_mri_fwhm(concMap, autoMask.data, options);
            obj.ttp = permute(ttpFit,[3,1,2]);

            % Take TTP as time from max of AIF function in seconds
            [~,aifmidx] = max(aifResult.conc(:));
            obj.ttp = obj.ttp - aifmidx;
            obj.ttp = obj.ttp*obj.tr/1000;

            if app.SVDButton.Value == 1
                obj.cbf = permute(cbfFit.svd.map,[3,1,2]);
                obj.mtt = permute(mttFit.svd,[3,1,2]);
            end

            if app.cSVDButton.Value == 1
                obj.cbf = permute(cbfFit.csvd.map,[3,1,2]);
                obj.mtt = permute(mttFit.csvd,[3,1,2]);
            end

            if app.oSVDButton.Value == 1
                obj.cbf = permute(cbfFit.osvd.map,[3,1,2]);
                obj.mtt = permute(mttFit.osvd,[3,1,2]);
            end

            obj.aif = aifResult.conc;
            obj.conc = concMap;
            obj.s0 = s0Map;
            obj.fwhm = fwhmFit;

            % Remove outliers
            obj.cbv(isnan(obj.cbv)) = 0;
            obj.cbv(isinf(obj.cbv)) = 0;
            obj.cbv(obj.cbv < 0) = 0;

            obj.cbvlc(isnan(obj.cbvlc)) = 0;
            obj.cbvlc(isinf(obj.cbvlc)) = 0;
            obj.cbvlc(obj.cbvlc < 0) = 0;

            obj.cbf(isnan(obj.cbf)) = 0;
            obj.cbf(isinf(obj.cbf)) = 0;
            obj.cbf(obj.cbf < 0) = 0;

            obj.mtt(isnan(obj.mtt)) = 0;
            obj.mtt(isinf(obj.mtt)) = 0;
            obj.mtt(obj.mtt < 0) = 0;

            obj.ttp(isnan(obj.ttp)) = 0;
            obj.ttp(isinf(obj.ttp)) = 0;
            obj.ttp(obj.ttp < 0) = 0;

            obj.k2map(isnan(obj.k2map)) = 0;
            obj.k2map(isinf(obj.k2map)) = 0;
            obj.k2map(obj.k2map < 0) = 0;

            % ROIs and points
            obj.ellipseROI = roiResult;
            app.roiPoints = [];
            app.roiPoints(:,1) = aifResult.voxels(:,1);
            app.roiPoints(:,2) = aifResult.voxels(:,2);
            app.roiPoints(:,3) = repmat(options.aif.nSlice,length(aifResult.voxels(:,1)),1);

            % Success
            obj.validFitFlag = true;

        end % DSCfit




        % ---------------------------------------------------------------------------------
        % Export the maps to GIFs
        % ---------------------------------------------------------------------------------
        function obj = exportGifDSC(obj, app, gifExportBase)

            % Scaling
            xScaling = obj.gifExportSize/obj.dimx;
            yScaling = obj.gifExportSize/obj.dimy/obj.aspectRatio;

            % Animated gif delay time
            delayTime = 2/obj.ns; % 2s

            % Create new directory
            ready = false;
            cnt = 1;
            while ~ready
                gifexportpath = strcat(gifExportBase,obj.tag,'DSC',filesep,num2str(cnt),filesep);
                if ~exist(gifexportpath, 'dir')
                    mkdir(gifexportpath);
                    ready = true;
                end
                cnt = cnt + 1;
            end

            % Write GIFs
            for slice = 1:obj.ns

                cbvim = uint8(round((255/app.CBVScaleEditField.Value)*obj.matrixInterpolate(squeeze(obj.cbv(slice,:,:).*obj.mask(slice,:,:)),[xScaling,yScaling],'makima')));
                cbfim = uint8(round((255/app.CBFScaleEditField.Value)*obj.matrixInterpolate(squeeze(obj.cbf(slice,:,:).*obj.mask(slice,:,:)),[xScaling,yScaling],'makima')));
                mttim = uint8(round((255/app.MTTScaleEditField.Value)*obj.matrixInterpolate(squeeze(obj.mtt(slice,:,:).*obj.mask(slice,:,:)),[xScaling,yScaling],'makima')));
                ttpim = uint8(round((255/app.TTPScaleEditField.Value)*obj.matrixInterpolate(squeeze(obj.ttp(slice,:,:).*obj.mask(slice,:,:)),[xScaling,yScaling],'makima')));
                cbvlcim = uint8(round((255/app.CBVlcScaleEditField.Value)*obj.matrixInterpolate(squeeze(obj.cbvlc(slice,:,:).*obj.mask(slice,:,:)),[xScaling,yScaling],'makima')));
                k2im = uint8(round((255/app.K2ScaleEditField.Value)*obj.matrixInterpolate(squeeze(obj.k2map(slice,:,:).*obj.mask(slice,:,:)),[xScaling,yScaling],'makima')));

                if slice == 1

                    imwrite(rot90(cbvim),app.cbvcmap,strcat(gifexportpath,filesep,'rCBV-',obj.tag,'.gif'),'DelayTime',delayTime,'LoopCount',inf);
                    imwrite(rot90(cbfim),app.cbfcmap,strcat(gifexportpath,filesep,'rCBF-',obj.tag,'.gif'),'DelayTime',delayTime,'LoopCount',inf);
                    imwrite(rot90(mttim),app.mttcmap,strcat(gifexportpath,filesep,'rMTT-',obj.tag,'.gif'),'DelayTime',delayTime,'LoopCount',inf);
                    imwrite(rot90(ttpim),app.ttpcmap,strcat(gifexportpath,filesep,'rTTP-',obj.tag,'.gif'),'DelayTime',delayTime,'LoopCount',inf);
                    imwrite(rot90(cbvlcim),app.cbvlccmap,strcat(gifexportpath,filesep,'rCBVlc-',obj.tag,'.gif'),'DelayTime',delayTime,'LoopCount',inf);
                    imwrite(rot90(k2im),app.k2cmap,strcat(gifexportpath,filesep,'K2-',obj.tag,'.gif'),'DelayTime',delayTime,'LoopCount',inf);

                else

                    imwrite(rot90(cbvim),app.cbvcmap,strcat(gifexportpath,filesep,'rCBV-',obj.tag,'.gif'),'WriteMode','append','DelayTime',delayTime);
                    imwrite(rot90(cbfim),app.cbfcmap,strcat(gifexportpath,filesep,'rCBF-',obj.tag,'.gif'),'WriteMode','append','DelayTime',delayTime);
                    imwrite(rot90(mttim),app.mttcmap,strcat(gifexportpath,filesep,'rMTT-',obj.tag,'.gif'),'WriteMode','append','DelayTime',delayTime);
                    imwrite(rot90(ttpim),app.ttpcmap,strcat(gifexportpath,filesep,'rTTP-',obj.tag,'.gif'),'WriteMode','append','DelayTime',delayTime);
                    imwrite(rot90(cbvlcim),app.cbvlccmap,strcat(gifexportpath,filesep,'rCBVlc-',obj.tag,'.gif'),'WriteMode','append','DelayTime',delayTime);
                    imwrite(rot90(k2im),app.k2cmap,strcat(gifexportpath,filesep,'K2-',obj.tag,'.gif'),'WriteMode','append','DelayTime',delayTime);

                end

            end

        end % exportGifDSC




        % ---------------------------------------------------------------------------------
        % Export the maps to DICOMs
        % ---------------------------------------------------------------------------------
        function obj = exportDCM(obj, dcmExportFolder)

            % Output maps
            dcmPar{1} = 'rCBV';
            dcmPar{2} = 'rCBF';
            dcmPar{3} = 'MTT';
            dcmPar{4} = 'TTP';
            dcmPar{5} = 'rCBVlc';
            dcmPar{6} = 'K2';

            dcmMap{1} = obj.cbv.*obj.mask;
            dcmMap{2} = obj.cbf.*obj.mask;
            dcmMap{3} = obj.mtt.*obj.mask;
            dcmMap{4} = obj.ttp.*obj.mask;
            dcmMap{5} = obj.cbvlc.*obj.mask;
            dcmMap{6} = obj.k2map.*obj.mask;

            dcmScale{1} = 1000;
            dcmScale{2} = 1000;
            dcmScale{3} = 1000;
            dcmScale{4} = 100;
            dcmScale{5} = 1000;
            dcmScale{6} = 1000;

            % Create new directories
            ready = false;
            cnt = 1;
            while ~ready
                outputFolderBase = strcat(dcmExportFolder,filesep,obj.tag,"DSC",filesep,num2str(cnt),filesep);
                if ~exist(outputFolderBase, 'dir')
                    mkdir(outputFolderBase);
                    ready = true;
                end
                cnt = cnt + 1;
            end

            for cnt = 1:length(dcmPar)
                outputFolder{cnt} = strcat(outputFolderBase,dcmPar{cnt});
                if ~exist(outputFolder{cnt}, 'dir')
                    mkdir(outputFolder{cnt});
                end
                delete(strcat(outputFolder{cnt},filesep,'*'));
            end
       
            % Export the maps
            for cnt = 1:length(dcmPar)

                seriesInstanceID = dicomuid;  % unique ID per parameter
                
                for slice = 1:size(dcmMap{cnt},1)

                    % Read the Dicom header
                    dcmHeader(slice) = obj.dcmInfo{slice}; %#ok<*AGROW>

                    % Changes some tags
                    dcmHeader(slice).ImageType = 'DERIVED\DSC\';
                    dcmHeader(slice).InstitutionName = 'Amsterdam UMC';
                    dcmHeader(slice).InstitutionAddress = 'Amsterdam, Netherlands';
                    dcmHeader(slice).StudyDescription = 'DSC mapping';

                    dcmHeader(slice).ProtocolName = dcmPar{cnt};
                    dcmHeader(slice).SequenceName = dcmPar{cnt};
                    dcmHeader(slice).SeriesDescription = 'DSC parameter map';
                    dcmHeader(slice).SeriesInstanceUID = seriesInstanceID;
                    dcmHeader(slice).RescaleSlope  = 1;
                    dcmHeader(slice).RescaleIntercept = 0;
                    dcmHeader(slice).NumberOfTemporalPositions = 1;
                    dcmHeader(slice).TemporalPositionIdentifier = 1;
                    
                    % Image position from sorted list
                    dcmHeader(slice).ImagePositionPatient = obj.position(obj.positionIndx(1,slice),:);
                    
                    fn = strcat('0000',num2str(slice));
                    fn = fn(size(fn,2)-4:size(fn,2));
                    fname = strcat(outputFolder{cnt},filesep,dcmPar{cnt},fn,'.dcm');
                    image = rot90(squeeze(cast(round(dcmScale{cnt}*dcmMap{cnt}(slice,:,:)),'uint16')));
                    dicomwrite(image, fname, dcmHeader(slice));

                end
            end

        end % exportDCM




    end % methods





    % ---------------------------------------------------------------------------------
    % Static methods
    % ---------------------------------------------------------------------------------
    methods (Static)



        % ---------------------------------------------------------------------------------
        % Resize an n-dimensional matrix
        % ---------------------------------------------------------------------------------
        function outputMatrix = matrixInterpolate(inputMatrix, scaling, varargin)

            % inputMatrx = n-dimensional matrix
            % scaling = scaling factor
            % varargin = 'linear', 'cubic' interpolation method

            N = ndims(inputMatrix);
            scaling = scaling(1:N);
            scaling(1,1:N) = scaling(:).';
            sz = size(inputMatrix);
            xvec = cell(1,N);
            yvec = cell(1,N);
            szy = nan(1,N);
            nonsing = true(1,N);

            for i = 1:N

                n = sz(i);

                if n==1 
                    nonsing(i) = 0;
                    szy(i) = 1;
                    continue
                end

                szy(i) = round(sz(i)*scaling(i));
                m = szy(i);

                xax = linspace(1/n/2, 1-1/n/2 ,n);
                xax = xax-.5;

                yax = linspace(1/m/2, 1-1/m/2 ,m);
                yax = yax-.5;

                xvec{i} = xax;
                yvec{i} = yax;

            end

            xvec = xvec(nonsing);
            yvec = yvec(nonsing);
            F = griddedInterpolant(xvec,squeeze(inputMatrix),varargin{:});
            outputMatrix = reshape(F(yvec),szy);

        end % matrixInterpolate




    end % Static methods




end