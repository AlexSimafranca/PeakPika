classdef Dataset
    properties (SetAccess = private)
        filenames  (:,1) {mustBeText}  = string.empty()
        data       (:,1) Pika.Spectrum = Pika.Spectrum.empty()
        background (1,1) Pika.Spectrum = Pika.Spectrum()
        peakFits   (:,1) Pika.PeakFit  = Pika.PeakFit.empty()
    end

    methods
        %%%---File Reading Functions---%%%
        function obj = readData(obj, path, options)
            arguments
                obj                     (1,1) Pika.Dataset
                path                    (1,:) {mustBeText}
                options.fileExt         (1,:) {mustBeText}    = [".txt", ".csv", ".dat", ".xls", ".xlsx", ".xlsm"]
                options.NamePrelimiter  (1,:) {mustBeText}    = ""
                options.NamePostlimiter (1,:) {mustBeText}    = ""
                options.HeaderLines     (1,1) {mustBeInteger} = 0
            end

            import Pika.Dir
            import Pika.Dataset
            import Pika.Spectrum

            % If path is specific file, read in data directly
            pathIsFile = contains(path, options.fileExt);
            if pathIsFile
                obj.filenames = Dir.getFilenameFromPath(path);
            else
                obj.filenames = Dir.read(path, options.fileExt);

                % Check that directory contained data files
                assert(~isempty(obj.filenames), ...
                       "[ Dataset > readData ] File read error: Could not read files from '" + path    + ...
                       "'. No files matching the specified filename parameters (Extensions: ['"        + ...
                       strjoin(options.fileExt, "', '") + "'], Prelimiter: '" + options.NamePrelimiter + ...
                       "', Postlimiter: '"  + options.NamePostlimiter + "') were found.")            
            end

            % Read in spectra from files
            numSpectra = length(obj.filenames);

            for i = 1:numSpectra
                % Construct full path to data file
                if pathIsFile
                    fullPath = path;
                else
                    fullPath = path + "/" + obj.filenames(i);
                end

                % Read in data
                newData = Spectrum.read(fullPath, "HeaderLines"    , options.HeaderLines     , ...
                                                  "NamePrelimiter" , options.NamePrelimiter  , ...
                                                  "NamePostlimiter", options.NamePostlimiter);

                % Append new data to data list
                obj = obj.appendData(newData);
            end
        end

        function obj = readBackground(obj, filePath, options)
            arguments
                obj                    (1,1) Pika.Dataset
                filePath               (1,:) {mustBeText}
                options.HeaderLines    (1,1) {mustBeInteger} = 0
            end

            import Pika.Spectrum

            % Check that file path has an extension
            [~, ~, ext] = fileparts(filePath);
            assert(~strcmp(ext, ""), ...
                   "[ Dataset > readBackground ] File error: Specified file path "   + ...
                   "does not have a file extension. Check that the path contains a " + ...
                   "specific filename, not just a directory.")

            % Read in background spectrum
            obj.background = Spectrum.read(filePath, "HeaderLines", options.HeaderLines);
        end
    
        %%%---Property Access Functions---%%%
        function obj = appendData(obj, spectrumArray)
            arguments
                obj           (1,1) Pika.Dataset
                spectrumArray (:,1) Pika.Spectrum
            end

            % Assign spectra unique IDs
            numSpectra = length(spectrumArray);
            excludeList = zeros(1, numSpectra);

            for i = 1:numSpectra
                newID = obj.assignSpectrumID("Exclude", excludeList(excludeList ~= 0));
                spectrumArray(i) = spectrumArray(i).setID(newID);
                
                % Update list of excluded IDs
                excludeList(i) = newID;
            end

            % Append spectra to end of data list
            obj.data = [obj.data; spectrumArray];
        end

        function ID = assignSpectrumID(obj, options)
            arguments
                obj             (1,1) Pika.Dataset
                options.Exclude (1,:) {mustBeInteger} = []
            end
            
            % Compile list of sample IDs for current data
            spectrumIDList = obj.querySpectrumIDs();

            % If spectrum ID list is empty, return early
            if isempty(spectrumIDList) && isempty(options.Exclude)
                ID = 1;
                return
            end

            % Find available ID to assign
            maxIndex = max([spectrumIDList, options.Exclude]) + 1;
            isAvailableID = ~ismember(1:maxIndex, spectrumIDList) & ~ismember(1:maxIndex, options.Exclude);

            ID = find(isAvailableID, 1);
        end

        function obj = setSpectrumName(obj, spectrumID, name)
            arguments
                obj        (1,1) Pika.Dataset
                spectrumID (1,1) {mustBeInteger, mustBePositive}
                name       (1,:) {mustBeText}
            end

            spectrumIndex = obj.findSpectrumIndexBySpectrumID(spectrumID);
            obj.data(spectrumIndex) = obj.data(spectrumIndex).setName(name);
        end

        %%%---Data Manipulation Functions---%%%
        function obj = subtractBackground(obj, scale, options)
            arguments
                obj             (1,1) Pika.Dataset
                scale           (1,1) {mustBeNumeric} = 1
                options.NameTag (1,:) {mustBeText}    = "\_bkg"
            end

            import Pika.Spectrum

            for i = 1:length(obj.data)
                % Calculate masks for common range
                [dataMask, backgroundMask] = Spectrum.findCommonRange(obj.data(i).X, obj.background.X);

                % Subtract background spectrum over common range
                newY = obj.data(i).Y(dataMask) - scale * obj.background.Y(backgroundMask);
                obj.data(i) = obj.data(i).updateY(newY, dataMask);

                % Append tag to end of spectrum name
                newName = obj.data(i).name + options.NameTag;
                obj.data(i) = obj.data(i).setName(newName);
            end
        end

        function obj = normalizeData(obj, maxRange, options)
            arguments 
                obj                  (1,1) Pika.Dataset
                maxRange             (:,2) {mustBeNumeric} = [NaN, NaN]
                options.MinRange     (:,2) {mustBeNumeric} = [NaN, NaN]
                options.SetMinToZero (1,:) {mustBeMember(options.SetMinToZero, ["Yes", "No"])} = "No"
            end

            % If only one max/min range specified, copy range for every spectrum
            numSpectra = length(obj.data);

            numMaxRanges = size(maxRange, 1);
            if numMaxRanges == 1
                maxRange = repmat(maxRange, numSpectra, 1);
            else 
                assert(numMaxRanges == numSpectra, ...
                       "[ Dataset > normalizeData ] Normalization error: Number of "    + ...
                       "bound value pairs for peak max range (" + string(numMaxRanges)  + ...
                       ") does not match the number of spectra in the dataset ("        + ...
                       string(numSpectra) + "). Either one range for all spectra or a " + ...
                       "range per spectrum is acceptable.")
            end

            numMinRanges = size(options.MinRange, 1);
            if numMinRanges == 1
                options.MinRange = repmat(options.MinRange, numSpectra, 1);
            else
                assert(numMinRanges == numSpectra, ...
                       "[ Dataset > normalizeData ] Normalization error: Number of "    + ...
                       "bound value pairs for range to zero (" + string(numMinRanges)   + ...
                       ") does not match the number of spectra in the dataset ("        + ...
                       string(numSpectra) + "). Either one range for all spectra or a " + ...
                       "range per spectrum is acceptable.")
            end

            % Normalize each spectrum
            for i = 1:numSpectra
                obj.data(i) = obj.data(i).normalize(maxRange(i,:), "MinRange"    , options.MinRange(i,:), ...
                                                                   "SetMinToZero", options.SetMinToZero);
            end
        end

        function obj = trimData(obj, dataRange, options)
            arguments
                obj               (1,1) Pika.Dataset
                dataRange         (1,2) {mustBeNumeric}
                options.RangeType (1,:) {mustBeMember(options.RangeType, ["Value", "Index"])} = "Value"
            end
            
            % Trim each spectrum in dataset
            numSpectra = length(obj.data);
            for i = 1:numSpectra
                obj.data(i) = obj.data(i).trimData(dataRange, "RangeType", options.RangeType);
            end

            % Trim background spectrum
            obj.background = obj.background.trimData(dataRange, "RangeType", options.RangeType);
        end

        function obj = maskData(obj, maskRange, spectrumID)
            arguments
                obj        (1,1) Pika.Dataset
                maskRange  (1,2) {mustBeNumeric}
                spectrumID (1,:) {mustBeInteger, mustBePositive} = obj.querySpectrumIDs()
            end

            for ID = spectrumID
                spectrumIndex = obj.findSpectrumIndexBySpectrumID(ID);
                obj.data(spectrumIndex) = obj.data(spectrumIndex).maskData(maskRange);
            end
        end

        % X-axis conversions
        function obj = convert2eV(obj)
            for i = 1:length(obj.data)
                obj.data(i) = obj.data(i).convert2eV();
            end
        end

        function obj = convert2nm(obj)
            for i = 1:length(obj.data)
                obj.data(i) = obj.data(i).convert2nm();
            end
        end

        %%%---Peak Fitting Functions---%%%
        function obj = initializePeakFits(obj, spectrumIDList)
            arguments
                obj            (1,1) Pika.Dataset
                spectrumIDList (1,:) {mustBeInteger, mustBePositive}
            end

            import Pika.PeakFit

            numFits = length(spectrumIDList);
            newPeakFits(numFits, 1) = PeakFit();

            for i = numFits:-1:1
                thisID = spectrumIDList(i);

                % If no corresponding PeakFit object, create one
                spectrumIndex = obj.findSpectrumIndexBySpectrumID(thisID);
                fitIndex = obj.findFitIndexBySpectrumID(thisID);

                if fitIndex == 0
                    thisSpectrum = obj.data(spectrumIndex);
                    newPeakFits(i) = PeakFit(thisSpectrum);
                    continue
                end

                % If PeakFit object already exists, remove preallocated memory
                newPeakFits(i) = [];
            end

            % Append new PeakFit objects to list of existing objects
            obj.peakFits = [obj.peakFits; newPeakFits];
        end

        function obj = addPeak(obj, peakCenter, sigma, options)
            arguments
                obj                    (1,1) Pika.Dataset
                peakCenter             (1,:) {mustBeNumeric}                                              = 0
                sigma                  (1,:) {mustBeNumeric, mustBePositive}                              = 1
                options.SpectrumIDList (1,:) {mustBeInteger, mustBePositive}                              = obj.querySpectrumIDs()
                options.PeakType       (1,:) {mustBeMember(options.PeakType, ["Gaussian", "Lorentzian"])} = "Gaussian"
                options.PeakHeight     (1,:) {mustBeNumeric}                                              = 1
            end

            import Pika.PeakFit

            % Check that all requested spectrum IDs exist
            obj.verifySpectrumIDs(options.SpectrumIDList);

            % If no existing PeakFit object for a spectrum, create one
            obj = obj.initializePeakFits(options.SpectrumIDList);

            % Check that valid number of each peak parameter provided
            numPeakCenter = length(peakCenter);

            numSigma = length(sigma);
            if numSigma == 1 && numPeakCenter ~= 1
                sigma = sigma * ones(1, numPeakCenter);
            else
                assert(numSigma == numPeakCenter, ...
                       "[ Dataset > addPeak ] Peak creation error: Number of sigma " + ...
                       "values (" + string(numSigma) + ") does not match number of " + ...
                       "peak centers (" + string(numPeakCenter) + ").")
            end

            numPeakHeight = length(options.PeakHeight);
            if numPeakHeight == 1 && numPeakCenter ~= 1
                options.PeakHeight = options.PeakHeight * ones(1, numPeakHeight);
            else
                assert(numPeakHeight == numPeakCenter, ...
                       "[ Dataset > addPeak ] Peak creation error: Number of peak height " + ...
                       "values (" + string(numPeakHeight) + ") does not match number of "  + ...
                       "peak centers (" + string(numPeakCenter) + ").")
            end


            % Add requested peak series to PeakFit object of each spectrum
            for ID = options.SpectrumIDList
                fitIndex = obj.findFitIndexBySpectrumID(ID);

                for i = 1:numPeakCenter
                    obj.peakFits(fitIndex) = obj.peakFits(fitIndex).addPeak(peakCenter(i), sigma(i),              ...
                                                                            "PeakType"  , options.PeakType      , ...
                                                                            "PeakHeight", options.PeakHeight(i));
                end
            end
        end

        function obj = addPeakSeries(obj, numPeaks, firstPeakCenter, peakSpacing, sigma, options)
            arguments
                obj                    (1,1) Pika.Dataset
                numPeaks               (1,1) {mustBeInteger, mustBeNonnegative}                           = 0
                firstPeakCenter        (1,1) {mustBeNumeric}                                              = 0
                peakSpacing            (1,1) {mustBeNumeric}                                              = 1
                sigma                  (1,1) {mustBeNumeric, mustBePositive}                              = 1
                options.SpectrumIDList (1,:) {mustBeInteger, mustBePositive}                              = obj.querySpectrumIDs()
                options.PeakType       (1,:) {mustBeMember(options.PeakType, ["Gaussian", "Lorentzian"])} = "Gaussian"
                options.PeakHeightList (1,:) {mustBeNumeric}                                              = ones(1, numPeaks)
            end
            
            import Pika.PeakFit

            % Check that all requested spectrum IDs exist
            obj.verifySpectrumIDs(options.SpectrumIDList);

            % If no existing PeakFit object for a spectrum, create one
            obj = obj.initializePeakFits(options.SpectrumIDList);

            % Add requested peak series to PeakFit object of each spectrum
            for ID = options.SpectrumIDList
                fitIndex = obj.findFitIndexBySpectrumID(ID);
                obj.peakFits(fitIndex) = obj.peakFits(fitIndex).addPeakSeries(numPeaks, firstPeakCenter, peakSpacing, sigma, ...
                                                                              "PeakType"      , options.PeakType       , ...
                                                                              "PeakHeightList", options.PeakHeightList);
            end
        end

        function obj = optimizeFit(obj, options)
            arguments
                obj
                options.FitIndex         (1,:) {mustBeInteger, mustBePositive}                         = 1:length(obj.peakFits)
                options.NumIterations    (1,1) {mustBeInteger, mustBePositive}                         = 100
                options.StepPeakCenter   (1,:) {mustBeMember(options.StepPeakCenter, ["Yes", "No"])}   = "Yes"
                options.PeakCenterStep   (1,1) {mustBeNumeric}                                         = 0.01
                options.StepPeakWidth    (1,:) {mustBeMember(options.StepPeakWidth, ["Yes", "No"])}    = "Yes"
                options.PeakWidthStep    (1,1) {mustBeNumeric}                                         = 0.01
                options.StepPeakHeight   (1,:) {mustBeMember(options.StepPeakHeight, ["Yes", "No"])}   = "Yes"
                options.PeakHeightStep   (1,1) {mustBeNumeric}                                         = 0.01
                options.SharePeakWidth   (1,1) {mustBeMember(options.SharePeakWidth, ["Yes", "No"])}   = "Yes"
                options.SharePeakSpacing (1,1) {mustBeMember(options.SharePeakSpacing, ["Yes", "No"])} = "Yes"
            end

            for i = options.FitIndex
                % Initialize fit progress plot
                X = NaN(options.NumIterations + 1, 1);
                X(1) = 0;
                Y = NaN(options.NumIterations + 1, 1);
                Y(1) = obj.peakFits(i).getFitRMSD();

                figure()
                h = plot(X, Y, "LineWidth", 2);

                % Plot settings
                set(gca, "FontSize", 15, "LineWidth", 1)
                title("PeakFit " + string(i))
                xlabel("Fit Iteration", "FontWeight", "bold")
                ylabel("RMSD", "FontWeight", "bold")
                xlim([0, options.NumIterations])
                ylim([0, inf])

                for j = 1:options.NumIterations
                    [obj.peakFits(i), RMSD] = obj.peakFits(i).stepPeaks("StepPeakCenter"  , options.StepPeakCenter   , ...
                                                                        "PeakCenterStep"  , options.PeakCenterStep   , ...
                                                                        "StepPeakWidth"   , options.StepPeakWidth    , ...
                                                                        "PeakWidthStep"   , options.PeakWidthStep    , ...
                                                                        "StepPeakHeight"  , options.StepPeakHeight   , ...
                                                                        "PeakHeightStep"  , options.PeakHeightStep   , ...
                                                                        "SharePeakWidth"  , options.SharePeakWidth   , ...
                                                                        "SharePeakSpacing", options.SharePeakSpacing);
                
                    % Update fit progress plot
                    h.XData(j+1) = j;
                    h.YData(j+1) = RMSD;

                    drawnow;
                end

                % Log fit parameters to command window
                disp("Peak Fit " + string(i) + ":")
                % disp(newline)
                disp("<strong>Fit RMSD: </strong>" + h.YData(end))
                disp("<strong>Peak Parameters: </strong>")
                disp(obj.peakFits(i).compileParameterTable())
            end
        end

        %%%---Plotting Functions---%%%
        function plotHandleList = plot(obj, axesHandle, options)
            arguments
                obj                       (1,1) Pika.Dataset
                axesHandle                (1,1) matlab.graphics.axis.Axes                             = matlab.graphics.axis.Axes()
                options.SpectrumIndexList (1,:) {mustBeInteger}                                       = 1:length(obj.data)
                options.NewFigure         (1,:) {mustBeMember(options.NewFigure, ["Yes", "No"])}      = "Yes"
                options.ShowBackground    (1,:) {mustBeMember(options.ShowBackground, ["Yes", "No"])} = "No"
            end

            import Pika.Spectrum

            % Plot requested spectra
            axesHandle = Spectrum.setAxes(axesHandle, "NewFigure", options.NewFigure);
            hold(axesHandle, "On")

            plotHandleList = gobjects(length(options.SpectrumIndexList), 1);
            for i = options.SpectrumIndexList
                plotHandleList(i) = obj.data(i).plot(axesHandle);
            end

            % If requested, plot background spectrum
            if strcmp(options.ShowBackground, "Yes")
                plotHandleList(end+1) = obj.background.plot(axesHandle);
            end
        end
    
        function plotHandles = plotFit(obj, axesHandle, options)
            arguments
                obj                    (1,1) Pika.Dataset
                axesHandle             (1,1) matlab.graphics.axis.Axes                        = matlab.graphics.axis.Axes()
                options.SpectrumIDList (1,:) {mustBeInteger, mustBePositive}                  = obj.queryFitSpectrumIDs()
                options.NewFigure      (1,:) {mustBeMember(options.NewFigure, ["Yes", "No"])} = "Yes"
                options.ShowPeaks      (1,:) {mustBeMember(options.ShowPeaks, ["Yes", "No"])} = "Yes"
            end

            import Pika.Spectrum

            % Plot fit objects
            plotHandles = cell(length(obj.peakFits), 1);
            for i = 1:length(options.SpectrumIDList)
                thisAxesHandle = Spectrum.setAxes(axesHandle, "NewFigure", options.NewFigure);

                thisID = options.SpectrumIDList(i);
                thisFitIndex = obj.findFitIndexBySpectrumID(thisID);
                plotHandles{i} = obj.peakFits(thisFitIndex).plot(thisAxesHandle, "ShowPeaks", options.ShowPeaks);
            end
        end

        %%%---Write To File Functions---%%%
        function writeToXLSX(obj, filename, options)
            arguments
                obj               (1,1) Pika.Dataset
                filename          (1,:) {mustBeText}
                options.RemoveNaN (1,1) {mustBeMember(options.RemoveNaN, ["Yes", "No"])} = "Yes"
            end
            
            for i = 1:length(obj.peakFits)
                obj.peakFits(i).writeToXLSX(filename, "PeakFitIndex", i, "RemoveNaN", options.RemoveNaN);
            end
        end

        %%%---Utility Functions---%%%
        function IDList = querySpectrumIDs(obj)
            numSpectra = length(obj.data);
            IDList = zeros(1, numSpectra);

            for i = 1:numSpectra
                IDList(i) = obj.data(i).ID;
            end
        end

        function fitIDList = queryFitSpectrumIDs(obj)
            numFits = length(obj.peakFits);
            fitIDList = zeros(numFits, 1);

            for i = 1:numFits
                fitIDList(i) = obj.peakFits(i).baseSpectrum.ID;
            end
        end

        function verifySpectrumIDs(obj, spectrumIDList)
            arguments
                obj            (1,1) Pika.Dataset
                spectrumIDList (1,:) {mustBeInteger, mustBePositive}
            end

            IDList = obj.querySpectrumIDs();

            % Check that number of requested IDs doesn't exceed available IDs
            assert(length(spectrumIDList) <= length(obj.data), ...
                   "[ Dataset > addPeakSeries ] Peak series addition error: Number " + ...
                   "of requested spectrum IDs exceeds number of available IDs ("     + ...
                   string(length(IDList)) + ").")

            % Check that all requested IDs correspond to spectra on object
            assert(all(ismember(spectrumIDList, IDList)), ...
                   "[ Dataset > addPeakSeries ] Peak series addition error: One or " + ...
                   "more of the requested spectrum IDs are invalid. Valid spectrum " + ...
                   "IDs are [" + join(string(IDList), ", ") + "].")
        end

        function index = findSpectrumIndexBySpectrumID(obj, spectrumID)
            arguments
                obj        (1,1) Pika.Dataset
                spectrumID (1,1) {mustBeInteger, mustBePositive}
            end

            index = 0;
            i = 1;
            while ~index && i <= length(obj.data)
                if obj.data(i).ID == spectrumID
                    index = i;
                    return
                end

                i = i + 1;
            end
        end

        function index = findFitIndexBySpectrumID(obj, spectrumID)
            arguments
                obj        (1,1) Pika.Dataset
                spectrumID (1,1) {mustBeInteger, mustBePositive}
            end

            index = 0;
            i = 1;
            while ~index && i <= length(obj.peakFits)
                if obj.peakFits(i).baseSpectrum.ID == spectrumID
                    index = i;
                    return
                end

                i = i + 1;
            end
        end
    end
end