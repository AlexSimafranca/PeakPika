% Peakfit.m

classdef PeakFit
    properties (SetAccess = private)
        baseSpectrum (1,1) Pika.Spectrum
        peakObjects  (:,1) {Pika.PeakFit.mustBeCellOfPeaks}
        numPeaks     (1,1) {mustBeInteger, mustBeNonnegative}
        fit          (:,1) Pika.Spectrum                    
    end

    methods (Static)
        function obj = PeakFit(dataSpectrum)
            arguments
                dataSpectrum (1,1) Pika.Spectrum = Pika.Spectrum()
            end

            % Store imported data spectrum as base spectrum
            obj.baseSpectrum = dataSpectrum;

            % Initialize empty peaks array
            obj.peakObjects = {};
            obj.numPeaks = 0;

            % Initialize empty fit spectrum
            obj.fit = dataSpectrum;
            obj.fit = obj.fit.setY([]);
            obj.fit = obj.fit.setName(obj.fit.name + "\_fit");
        end

        %%%---Utility Functions---%%%
        function mustBeCellOfPeaks(argument)
            arguments
                argument (:,1) cell
            end
            
            if isempty(argument)
                return
            end
            
            for i = 1:length(argument)
                if ~(isa(argument{i}, "Pika.Peak") || isa(argument{i}, "Pika.PeakSeries"))
                    error("All elements in the cell array must be of type 'Peak' or 'PeakSeries'.");
                end
            end
        end

        %%%---Helper Functions---%%%
        function totalYDiffStruct = calcTotalYDiffStruct(YDiffStructList)
            arguments
                YDiffStructList (:,1) struct
            end
            
            % Sum all Y difference vectors
            parameterKeys = ["peakCenter", "peakWidth", "peakHeight", "peakSpacing"];
            numStructs = length(YDiffStructList);

            totalYDiffStruct = YDiffStructList(1);
            for i = 2:numStructs
                for j = 1:length(parameterKeys)
                    totalYDiffStruct.(parameterKeys(j)) = totalYDiffStruct.(parameterKeys(j)) + YDiffStructList(i).(parameterKeys(j));
                end
            end
        end
    end

    methods
        %%%---Data Manipulation Functions---%%%
        % X-axis conversions
        function obj = convert2eV(obj)
            obj.baseSpectrum = obj.baseSpectrum.convert2eV();
            obj.fit = obj.fit.convert2eV();

            for i = 1:length(obj.peakObjects)
                obj.peakObjects{i} = obj.peakObjects{i}.convert2eV();
            end
        end

        function obj = convert2nm(obj)
            obj.baseSpectrum = obj.baseSpectrum.convert2nm();
            obj.fit = obj.fit.convert2nm();
            
            for i = 1:length(obj.peakObjects)
                obj.peakObjects{i} = obj.peakObjects{i}.convert2nm();
            end
        end

        %%%---Peak Fitting Functions---%%%
        function obj = addPeak(obj, peakCenter, sigma, options)
            arguments
                obj                    (1,1) Pika.PeakFit
                peakCenter             (1,1) {mustBeNumeric}                                              = 0
                sigma                  (1,1) {mustBeNumeric, mustBePositive}                              = 1
                options.PeakType       (1,:) {mustBeMember(options.PeakType, ["Gaussian", "Lorentzian"])} = "Gaussian"
                options.PeakHeight     (1,:) {mustBeNumeric}                                              = 1
            end

            import Pika.Peak

            % Generate new peak
            X = obj.baseSpectrum.X;
            newPeak = Peak(X, peakCenter, sigma, "PeakType"  , options.PeakType   , ...
                                                 "PeakHeight", options.PeakHeight);

            % Append new peak to peaks array
            obj.peakObjects(end+1, 1) = {newPeak};

            % Calculate total fit
            obj = obj.calculateFit();

            % Update number of peaks
            obj.numPeaks = obj.numPeaks + 1;
        end

        function obj = addPeakSeries(obj, numPeaks, firstPeakCenter, peakSpacing, sigma, options)
            arguments
                obj                    (1,1) Pika.PeakFit
                numPeaks               (1,1) {mustBeInteger, mustBeNonnegative}                           = 0
                firstPeakCenter        (1,1) {mustBeNumeric}                                              = 0
                peakSpacing            (1,1) {mustBeNumeric}                                              = 1
                sigma                  (1,1) {mustBeNumeric, mustBePositive}                              = 1
                options.PeakType       (1,:) {mustBeMember(options.PeakType, ["Gaussian", "Lorentzian"])} = "Gaussian"
                options.PeakHeightList (1,:) {mustBeNumeric}                                              = ones(1, numPeaks)
            end

            import Pika.PeakSeries

            % Generate new peak series
            X = obj.baseSpectrum.X;
            newPeakSeries = PeakSeries(X, numPeaks, firstPeakCenter, peakSpacing, sigma, ...
                                       "PeakType"      , options.PeakType       , ...
                                       "PeakHeightList", options.PeakHeightList);

            % Append new peak series to peaks array
            obj.peakObjects(end+1, 1) = {newPeakSeries};

            % Calculate total fit
            obj = obj.calculateFit();

            % Update number of peaks
            obj.numPeaks = obj.numPeaks + newPeakSeries.numPeaks;
        end

        function obj = calculateFit(obj)
            import Pika.Spectrum

            if isempty(obj.peakObjects)
                obj.fit = Spectrum();
                return
            end

            % Compile Y data from each peak object
            Y = zeros(length(obj.baseSpectrum.X), 1);
            for i = 1:length(obj.peakObjects)
                Y = Y + obj.peakObjects{i}.calculateFit();
            end

            % Assign fit object
            obj.fit = obj.baseSpectrum;
            obj.fit = obj.fit.setY(Y);
            obj.fit = obj.fit.setName(obj.fit.name + "\_fit");
        end

        function [obj, RMSD] = stepPeaks(obj, options)
            arguments
                obj                      (1,1) Pika.PeakFit
                options.StepPeakCenter   (1,:) {mustBeMember(options.StepPeakCenter, ["Yes", "No"])}   = "Yes"
                options.PeakCenterStep   (1,1) {mustBeNumeric}                                         = 1
                options.StepPeakWidth    (1,:) {mustBeMember(options.StepPeakWidth, ["Yes", "No"])}    = "Yes"
                options.PeakWidthStep    (1,1) {mustBeNumeric}                                         = 0.01
                options.StepPeakHeight   (1,:) {mustBeMember(options.StepPeakHeight, ["Yes", "No"])}   = "Yes"
                options.PeakHeightStep   (1,1) {mustBeNumeric}                                         = 0.01
                options.SharePeakWidth   (1,1) {mustBeMember(options.SharePeakWidth, ["Yes", "No"])}   = "Yes"
                options.SharePeakSpacing (1,1) {mustBeMember(options.SharePeakSpacing, ["Yes", "No"])} = "Yes"
            end

            import Pika.Stats

            % Step each peak to increase agreement with base spectrum
            numPeakObjects = length(obj.peakObjects);
            for i = 1:numPeakObjects
                switch class(obj.peakObjects{i})
                    case "Pika.Peak"
                        obj = stepSinglePeak(obj, i, "StepPeakCenter", options.StepPeakCenter , ...
                                                     "PeakCenterStep", options.PeakCenterStep , ...
                                                     "StepPeakWidth" , options.StepPeakWidth  , ...
                                                     "PeakWidthStep" , options.PeakWidthStep  , ...
                                                     "StepPeakHeight", options.StepPeakHeight , ...
                                                     "PeakHeightStep", options.PeakHeightStep);
                    case "Pika.PeakSeries"
                        if strcmp(options.SharePeakWidth, "No") && strcmp(options.SharePeakSpacing, "No")
                            obj = stepPeakSeries(obj, i, "StepPeakCenter", options.StepPeakCenter , ...
                                                         "PeakCenterStep", options.PeakCenterStep , ...
                                                         "StepPeakWidth" , options.StepPeakWidth  , ...
                                                         "PeakWidthStep" , options.PeakWidthStep  , ...
                                                         "StepPeakHeight", options.StepPeakHeight , ...
                                                         "PeakHeightStep", options.PeakHeightStep);
                        else
                            obj = stepSharedPeakSeries(obj, i, "StepPeakCenter"  , options.StepPeakCenter   , ...
                                                               "PeakCenterStep"  , options.PeakCenterStep   , ...
                                                               "StepPeakWidth"   , options.StepPeakWidth    , ...
                                                               "PeakWidthStep"   , options.PeakWidthStep    , ...
                                                               "StepPeakHeight"  , options.StepPeakHeight   , ...
                                                               "PeakHeightStep"  , options.PeakHeightStep   , ...
                                                               "SharePeakWidth"  , options.SharePeakWidth   , ...
                                                               "SharePeakSpacing", options.SharePeakSpacing);
                        end
                end
            end

            % Calculate root-mean-square difference between data and current fit
            RMSD = Stats.calculateRMSD(obj.baseSpectrum.Y, obj.fit.Y);
        end

        function obj = stepSinglePeak(obj, peakIndex, options)
            arguments
                obj                      (1,1) Pika.PeakFit
                peakIndex                (1,1) {mustBeInteger, mustBePositive}
                options.StepPeakCenter   (1,:) {mustBeMember(options.StepPeakCenter, ["Yes", "No"])}   = "Yes"
                options.PeakCenterStep   (1,1) {mustBeNumeric}                                         = 1
                options.StepPeakWidth    (1,:) {mustBeMember(options.StepPeakWidth, ["Yes", "No"])}    = "Yes"
                options.PeakWidthStep    (1,1) {mustBeNumeric}                                         = 0.01
                options.StepPeakHeight   (1,:) {mustBeMember(options.StepPeakHeight, ["Yes", "No"])}   = "Yes"
                options.PeakHeightStep   (1,1) {mustBeNumeric}                                         = 0.01
            end

            % Step requested peak parameters
            YDiffStruct = obj.peakObjects{peakIndex}.stepPeakParameters("StepPeakCenter", options.StepPeakCenter , ...
                                                                        "PeakCenterStep", options.PeakCenterStep , ...
                                                                        "StepPeakWidth" , options.StepPeakWidth  , ...
                                                                        "PeakWidthStep" , options.PeakWidthStep  , ...
                                                                        "StepPeakHeight", options.StepPeakHeight , ...
                                                                        "PeakHeightStep", options.PeakHeightStep);

            % Calculate updated fit Y data
            fitYStruct = obj.calcFitYStruct(YDiffStruct);

            % Evaluate possible parameter steps
            [parameterToStep, stepIncrement] = obj.decideParameterToStep(fitYStruct, ...
                                                                         "PeakCenterStep", options.PeakCenterStep , ...
                                                                         "PeakWidthStep" , options.PeakWidthStep  , ...
                                                                         "PeakHeightStep", options.PeakHeightStep);

            if stepIncrement ~= 0
                % Step peak
                obj.peakObjects{peakIndex} = obj.peakObjects{peakIndex}.stepPeak(parameterToStep, stepIncrement);

                % Recalculate fit
                obj = obj.calculateFit();
            end
        end

        function obj = stepPeakSeries(obj, peakIndex, options)
            arguments
                obj                    (1,1) Pika.PeakFit
                peakIndex              (1,1) {mustBeInteger, mustBePositive}
                options.StepPeakCenter (1,:) {mustBeMember(options.StepPeakCenter, ["Yes", "No"])} = "Yes"
                options.PeakCenterStep (1,1) {mustBeNumeric}                                       = 1
                options.StepPeakWidth  (1,:) {mustBeMember(options.StepPeakWidth, ["Yes", "No"])}  = "Yes"
                options.PeakWidthStep  (1,1) {mustBeNumeric}                                       = 0.01
                options.StepPeakHeight (1,:) {mustBeMember(options.StepPeakHeight, ["Yes", "No"])} = "Yes"
                options.PeakHeightStep (1,1) {mustBeNumeric}                                       = 0.01
            end

            % Step requested peak parameters
            YDiffStructList = obj.peakObjects{peakIndex}.stepPeakParameters("StepPeakCenter", options.StepPeakCenter , ...
                                                                            "PeakCenterStep", options.PeakCenterStep , ...
                                                                            "StepPeakWidth" , options.StepPeakWidth  , ...
                                                                            "PeakWidthStep" , options.PeakWidthStep  , ...
                                                                            "StepPeakHeight", options.StepPeakHeight , ...
                                                                            "PeakHeightStep", options.PeakHeightStep);

            % Step each peak in direction of better agreement
            for i = 1:length(YDiffStructList)
                % Calculate updated fit Y data
                fitYStruct = obj.calcFitYStruct(YDiffStructList(i));

                % Evaluate possible parameter steps
                [parameterToStep, stepIncrement] = obj.decideParameterToStep(fitYStruct, ...
                                                                             "PeakCenterStep", options.PeakCenterStep , ...
                                                                             "PeakWidthStep" , options.PeakWidthStep  , ...
                                                                             "PeakHeightStep", options.PeakHeightStep);
                
                if stepIncrement ~= 0
                    % Step peak
                    obj.peakObjects{peakIndex} = obj.peakObjects{peakIndex}.stepPeak(parameterToStep, stepIncrement, i);

                    % Recalculate fit
                    obj = obj.calculateFit();
                end
            end        
        end

        function obj = stepSharedPeakSeries(obj, peakIndex, options)
            arguments
                obj                      (1,1) Pika.PeakFit
                peakIndex                (1,1) {mustBeInteger, mustBePositive}
                options.StepPeakCenter   (1,:) {mustBeMember(options.StepPeakCenter, ["Yes", "No"])}   = "Yes"
                options.PeakCenterStep   (1,1) {mustBeNumeric}                                         = 1
                options.StepPeakWidth    (1,:) {mustBeMember(options.StepPeakWidth, ["Yes", "No"])}    = "Yes"
                options.PeakWidthStep    (1,1) {mustBeNumeric}                                         = 0.01
                options.StepPeakHeight   (1,:) {mustBeMember(options.StepPeakHeight, ["Yes", "No"])}   = "Yes"
                options.PeakHeightStep   (1,1) {mustBeNumeric}                                         = 0.01
                options.StepPeakSpacing  (1,:) {mustBeMember(options.StepPeakSpacing, ["Yes", "No"])}  = "Yes"
                options.PeakSpacingStep  (1,1) {mustBeNumeric}                                         = 0.01
                options.SharePeakWidth   (1,1) {mustBeMember(options.SharePeakWidth, ["Yes", "No"])}   = "Yes"
                options.SharePeakSpacing (1,1) {mustBeMember(options.SharePeakSpacing, ["Yes", "No"])} = "Yes"
            end
            
            import Pika.PeakFit

            %---Phase 1 - Step Peaks Together---%
            % Step requested peak parameters
            YDiffStructList = obj.peakObjects{peakIndex}.stepSharedPeakParameters("StepPeakCenter"  , options.StepPeakCenter   , ...
                                                                                  "PeakCenterStep"  , options.PeakCenterStep   , ...
                                                                                  "StepPeakWidth"   , options.StepPeakWidth    , ...
                                                                                  "PeakWidthStep"   , options.PeakWidthStep    , ...
                                                                                  "StepPeakHeight"  , options.StepPeakHeight   , ...
                                                                                  "PeakHeightStep"  , options.PeakHeightStep   , ...
                                                                                  "StepPeakSpacing" , options.StepPeakSpacing  , ...
                                                                                  "PeakSpacingStep" , options.PeakSpacingStep  , ...
                                                                                  "SharePeakWidth"  , options.SharePeakWidth   , ...
                                                                                  "SharePeakSpacing", options.SharePeakSpacing);
            
            % Calculate average of Y difference vectors
            totalYDiffStruct = PeakFit.calcTotalYDiffStruct(YDiffStructList);
            
            % Calculate fit Y vectors for peak series
            totalFitYStruct = calcFitYStruct(obj, totalYDiffStruct);

            % Evaluate possible parameter steps
            [parameterToStep, stepIncrement] = obj.decideParameterToStep(totalFitYStruct, ...
                                                                         "PeakCenterStep" , options.PeakCenterStep  , ...
                                                                         "PeakWidthStep"  , options.PeakWidthStep   , ...
                                                                         "PeakHeightStep" , options.PeakHeightStep  , ...
                                                                         "PeakSpacingStep", options.PeakSpacingStep);
            
            if stepIncrement ~= 0
                % Step peak
                obj.peakObjects{peakIndex} = obj.peakObjects{peakIndex}.stepSharedPeaks(parameterToStep, stepIncrement);

                % Recalculate fit
                obj = obj.calculateFit();
            end

            %---Phase 2 - Step Peaks Individually---%
            % Translate input parameters
            if strcmp(options.SharePeakWidth, "Yes")
                options.StepPeakWidth = "No";
            end

            if strcmp(options.SharePeakSpacing, "Yes")
                options.StepPeakCenter = "No";
            end

            % Step peak series as individual peaks
            obj = stepPeakSeries(obj, peakIndex, "StepPeakCenter", options.StepPeakCenter , ...
                                                 "PeakCenterStep", options.PeakCenterStep , ...
                                                 "StepPeakWidth" , options.StepPeakWidth  , ...
                                                 "PeakWidthStep" , options.PeakWidthStep  , ...
                                                 "StepPeakHeight", options.StepPeakHeight , ...
                                                 "PeakHeightStep", options.PeakHeightStep);
        end

        function fitYStruct = calcFitYStruct(obj, YDiffStruct)
            arguments
                obj         (1,1) Pika.PeakFit
                YDiffStruct (1,1) struct
            end

            parameterKeys = fieldnames(YDiffStruct);
            fitYStruct = YDiffStruct;
            
            for i = 1:3
                if ~isempty(fitYStruct.(parameterKeys{i}))
                    fitYStruct.(parameterKeys{i}) = fitYStruct.(parameterKeys{i}) + obj.fit.Y * [1,1];
                end
            end
        end

        function [parameterToStep, stepIncrement] = decideParameterToStep(obj, YStruct, options)
            arguments
                obj                     (1,1) Pika.PeakFit
                YStruct                 (1,1) struct
                options.PeakCenterStep  (1,1) {mustBeNumeric} = 0
                options.PeakWidthStep   (1,1) {mustBeNumeric} = 0     
                options.PeakHeightStep  (1,1) {mustBeNumeric} = 0
                options.PeakSpacingStep (1,1) {mustBeNumeric} = 0
            end

            import Pika.Stats

            % Allocate storage for root-mean-square difference (RMSD) values
            RMSDSet = NaN(1,9);

            % Store current RMSD as last value
            RMSDSet(end) = Stats.calculateRMSD(obj.baseSpectrum.Y, obj.fit.Y);

            % Scan parameters for change to fit
            parameterKeys = ["peakCenter", "peakWidth", "peakHeight", "peakSpacing"];
            structParameterKeys = fieldnames(YStruct);

            for i = 1:length(parameterKeys)
                if any(ismember(structParameterKeys, parameterKeys(i)))
                    % Retrieve stepped Y data
                    thisYPair = YStruct.(parameterKeys(i));
                    if isempty(thisYPair)
                        continue
                    end

                    % Calculate RMSD with base spectrum
                    RMSDSet(2*(i-1) + 1) = Stats.calculateRMSD(obj.baseSpectrum.Y, thisYPair(:,1));
                    RMSDSet(2*i) = Stats.calculateRMSD(obj.baseSpectrum.Y, thisYPair(:,2));
                end
            end

            % Select step that minimizes RMSD
            parameterIndex = find(RMSDSet == min(RMSDSet), 1);
            parameterToStep = "";
            stepIncrement = 0;

            if parameterIndex == 7 
                return
            end

            switch ceil(parameterIndex / 2)
                case 1
                    parameterToStep = "peakCenter";
                    stepIncrement = options.PeakCenterStep;
                case 2
                    parameterToStep = "peakWidth";
                    stepIncrement = options.PeakWidthStep;
                case 3
                    parameterToStep = "peakHeight";
                    stepIncrement = options.PeakHeightStep;
                case 4
                    parameterToStep = "peakSpacing";
                    stepIncrement = options.PeakSpacingStep;
            end
    
            if mod(parameterIndex, 2) == 0
                stepIncrement = -stepIncrement;
            end
        end

        %%%---Plotting Functions---%%%
        function plotHandles = plot(obj, axesHandle, options)
            arguments
                obj               (1,1) Pika.PeakFit
                axesHandle        (1,1) matlab.graphics.axis.Axes                        = matlab.graphics.axis.Axes()
                options.NewFigure (1,:) {mustBeMember(options.NewFigure, ["Yes", "No"])} = "Yes"
                options.ShowPeaks (1,:) {mustBeMember(options.ShowPeaks, ["Yes", "No"])} = "Yes"
            end
            
            import Pika.Spectrum

            axesHandle = Spectrum.setAxes(axesHandle, "NewFigure", options.NewFigure);

            % Plot base spectrum
            plotHandles = struct();
            plotHandles.base = obj.baseSpectrum.plot(axesHandle);

            % Plot fit
            plotHandles.fit = obj.fit.plot(axesHandle);

            % Plot peaks, if requested
            if strcmp(options.ShowPeaks, "Yes")
                peakHandleList = gobjects(obj.numPeaks, 1);

                j = 1;
                for i = 1:length(obj.peakObjects)
                    numNewHandles = obj.peakObjects{i}.countPeaks();
                    jRange = [j, j + (numNewHandles - 1)];
                    peakHandleList(jRange(1):jRange(2)) = obj.peakObjects{i}.plot(axesHandle);
                    
                    % Set color of peak/peak series
                    for k = jRange(1):jRange(2)
                        set(peakHandleList(k), "Color", Spectrum.getColor(i))
                    end

                    % Label peak w/ peak number
                    [labelX, labelY] = calcLabelPosition(obj, i);
                    text(labelX, labelY, string(i), "FontSize", 15, ...
                                                    "HorizontalAlignment", "center")

                    j = j + numNewHandles;
                end        
            else 
                peakHandleList = [];
            end
            plotHandles.peaks = peakHandleList;

            % Plot settings
            legend([obj.baseSpectrum.name, obj.fit.name])
            set(plotHandles.base, "Color", "k")
            set(plotHandles.fit, "Color", "r", "LineStyle", "--")
        end

        %%%---Write To File Functions---%%%
        function writeToXLSX(obj, filename, options)
            arguments
                obj                  (1,1) Pika.PeakFit
                filename             (1,:) {mustBeText}
                options.Sheet        (1,:) {mustBeText}                                     = ""
                options.PeakFitIndex (1,1) {mustBeInteger, mustBePositive}                  = 1
                options.RemoveNaN    (1,1) {mustBeMember(options.RemoveNaN, ["Yes", "No"])} = "Yes"
            end

            % Set sheet name
            if strcmp(options.Sheet, "")
                options.Sheet = "peakfit" + options.PeakFitIndex;
            end

            % Peak fit RMSD
            writecell({"RMSD:", obj.getFitRMSD()}, filename, "Sheet", options.Sheet)
            
            % Peak parameters
            peakParameters = obj.compileParameterTable();

            startRow = 2;
            endRow = startRow + size(peakParameters, 1);
            rowRange = "A" + string(startRow);

            writetable(peakParameters, filename, "Sheet", options.Sheet, "Range", rowRange);
            
            % Peak Spectra
            baseSpectrumTable = obj.baseSpectrum.getSpectrumData("ReturnType", "Table", ...
                                                                 "RemoveNaN", options.RemoveNaN);
            baseSpectrumTable.Properties.VariableNames = ["Base Spectrum - X", "Base Spectrum - Y"];

            fitTable = obj.fit.getSpectrumData("ReturnType", "Table", ...
                                               "RemoveNaN", options.RemoveNaN);
            fitTable.Properties.VariableNames = ["Total Fit - X", "Total Fit - Y"];

            peakSpectra = obj.compileSpectrumData("RemoveNaN", options.RemoveNaN);

            startRow = endRow + 2;
            rowRange = "A" + string(startRow);

            writetable([baseSpectrumTable, fitTable, peakSpectra], filename, "Sheet", options.Sheet, "Range", rowRange);
        end

        %%%---Utility Functions---%%%
        function bool = isEmpty(obj)
            bool = obj.baseSpectrum.isEmpty();
        end

        function [labelX, labelY] = calcLabelPosition(obj, peakObjectIndex)
            arguments
                obj             (1,1) Pika.PeakFit
                peakObjectIndex (1,1) {mustBeInteger, mustBePositive}
            end

            import Pika.Spectrum

            % Retrieve peak position and height
            peakObject = obj.peakObjects{peakObjectIndex};
            switch class(peakObject)
                case "Pika.Peak"
                    peakPosition = peakObject.peakCenter;
                    peakY = peakObject.peakSpectrum.Y;
                case "Pika.PeakSeries"
                    peakPosition = peakObject.firstPeakCenter;
                    peakY = peakObject.peaks(1).peakSpectrum.Y;
                otherwise
                    disp("[ PeakFit > calcLabelPosition ] Peak label error: Peak object " + ...
                         string(peakObjectIndex) + " is not of an acceptable type. "      + ...
                         "Acceptable types are 'Pika.Peak' and Pika.PeakSeries'.")
                    labelX = 0;
                    labelY = 0;
                    return
            end

            % Bound the peak label position by the edges of the spectrum
            XBounds = [min(obj.baseSpectrum.X), max(obj.baseSpectrum.X)];
            YBounds = [min(obj.baseSpectrum.Y), max(obj.baseSpectrum.Y)];
            XRange = XBounds(2) - XBounds(1);
            YRange = YBounds(2) - YBounds(1);

            if peakPosition < XBounds(1)
                labelX = XBounds(1) + 0.02 * XRange;
                labelY = 1.05 * peakY(Spectrum.findXValue(obj.baseSpectrum.X, XBounds(1)));
            elseif peakPosition > XBounds(2)
                labelX = XBounds(2) + 0.02 * XRange;
                labelY = 1.05 * peakY(Spectrum.findXValue(obj.baseSpectrum.X, labelX));
            else
                labelX = peakPosition;
                labelY = max(peakY) + 0.05 * YRange;
            end
        end

        function RMSD = getFitRMSD(obj)
            import Pika.Stats

            RMSD = Stats.calculateRMSD(obj.baseSpectrum.Y, obj.fit.Y);
        end

        %%%---Helper Functions---%%%
        function parameterTable = compileParameterTable(obj)
            parameterTable = table();

            rowIndex = 1;
            for i = 1:length(obj.peakObjects)

                switch class(obj.peakObjects{i})
                    case "Pika.Peak"
                        parameterTable(rowIndex, :) = [{i}, {1}, obj.peakObjects{i}.getPeakParameters()];
                        
                        rowIndex = rowIndex + 1;

                    case "Pika.PeakSeries"
                        peakSeriesParameters = obj.peakObjects{i}.getPeakParameters();

                        for j = 1:size(peakSeriesParameters, 1)
                            parameterTable(rowIndex, :) = [{i}, peakSeriesParameters(j,:)];
                            rowIndex = rowIndex + 1;
                        end
                end
            end

            parameterTable.Properties.VariableNames = ["Peak Object #"     , ...
                                                       "Series Peak #"     , ...
                                                       "Peak Type"         , ...
                                                       "Peak Center"       , ...
                                                       "Peak Width (Sigma)", ...
                                                       "Peak Height"      ];        
        end

        function spectrumTable = compileSpectrumData(obj, options)
            arguments
                obj               (1,1) Pika.PeakFit
                options.RemoveNaN (1,1) {mustBeMember(options.RemoveNaN, ["Yes", "No"])} = "Yes"
            end

            spectrumTable = table();

            columnIndexPair = [1, 2];
            columnNames = strings(1, 2 * obj.numPeaks);
            for i = 1:length(obj.peakObjects)
                spectrumData = obj.peakObjects{i}.getSpectrumData("RemoveNaN", options.RemoveNaN);
                thisNumPeaks = size(spectrumData, 2) / 2;

                for j = 1:thisNumPeaks
                    spectrumTable(:, columnIndexPair(1)) = table(spectrumData(:, 2*(j-1) + 1));
                    spectrumTable(:, columnIndexPair(2)) = table(spectrumData(:, 2*j));
                    
                    columnNames(columnIndexPair) = ...
                        {"Peak " + string(i) + "." + string(j) + " - X", ...
                         "Peak " + string(i) + "." + string(j) + " - Y"};

                    columnIndexPair = columnIndexPair + 2;
                end
            end

            spectrumTable.Properties.VariableNames = columnNames;
        end    
    end
end 

