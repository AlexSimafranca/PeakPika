classdef PeakSeries
    properties
        peaks           (:,1) Pika.Peak
        numPeaks        (1,1) {mustBeInteger, mustBeNonnegative}
        firstPeakCenter (1,1) {mustBeNumeric}
        peakSpacing     (1,1) {mustBeNumeric}
        sigma           (1,1) {mustBeNumeric, mustBeNonnegative}
        peakHeightList  (1,:) {mustBeNumeric}
        peakType        (1,:) {mustBeMember(peakType, ["Gaussian", "Lorentzian"])} = "Gaussian"
    end

    methods (Static)
        function obj = PeakSeries(X, numPeaks, firstPeakCenter, peakSpacing, sigma, options)
            arguments
                X                      (:,1) {mustBeNumeric}                                              = []
                numPeaks               (1,1) {mustBeInteger, mustBeNonnegative}                           = 0
                firstPeakCenter        (1,1) {mustBeNumeric}                                              = 0
                peakSpacing            (1,1) {mustBeNumeric}                                              = 10
                sigma                  (1,1) {mustBeNumeric, mustBeNonnegative}                           = 1
                options.PeakHeightList (1,:) {mustBeNumeric}                                              = ones(1, numPeaks)
                options.PeakType       (1,:) {mustBeMember(options.PeakType, ["Gaussian", "Lorentzian"])} = "Gaussian"
            end
                
            import Pika.Peak
            import Pika.PeakShape

            obj.peakType = options.PeakType;
            
            % Store peak parameters on object
            obj.numPeaks = numPeaks;
            obj.firstPeakCenter = firstPeakCenter;
            obj.peakSpacing = peakSpacing;
            obj.sigma = sigma;
            obj.peakHeightList = options.PeakHeightList;
            obj.peakType = options.PeakType;

            % Check that peak heights are provided for every peak
            assert(length(options.PeakHeightList) == numPeaks, ...
                   "[ PeakSeries > PeakSeries] Peak series constructor error: Number of " + ...
                   "peak heights provide does not match the number of requested peaks.")

            % Create set of equally spaced peaks
            if isempty(X)
                obj.peaks = Pika.Peak.empty();
                return
            end

            obj.peaks(numPeaks, 1) = Peak();
            for i = 1:numPeaks
                thisPeakCenter = firstPeakCenter + (i-1) * peakSpacing;
                obj.peaks(i) = Peak(X, thisPeakCenter, sigma, "PeakType"  , options.PeakType          , ...
                                                              "PeakHeight", options.PeakHeightList(i));
            end
        end
    end

    methods
        %%%---Data Manipulation Functions---%%%
        % X-axis conversions
        function obj = convert2eV(obj)
            for i = 1:length(obj.peaks)
                obj.peaks(i) = obj.peaks(i).convert2eV();
            end
        end

        function obj = convert2nm(obj)
            for i = 1:length(obj.peaks)
                obj.peaks(i) = obj.peaks(i).convert2nm();
            end
        end
    
        %%%---Plotting Functions---%%%
        function plotHandleList = plot(obj, axesHandle, options)
            arguments
                obj               (1,1) Pika.PeakSeries
                axesHandle        (1,1) matlab.graphics.axis.Axes                        = matlab.graphics.axis.Axes()
                options.NewFigure (1,:) {mustBeMember(options.NewFigure, ["Yes", "No"])} = "Yes"
            end

            import Pika.Spectrum

            axesHandle = Spectrum.setAxes(axesHandle, "NewFigure", options.NewFigure);

            plotHandleList = gobjects(obj.numPeaks, 1);
            for i = 1:obj.numPeaks
                plotHandleList(i) = obj.peaks(i).plot(axesHandle);
            end
        end

        %%%---Peak Fitting Functions---%%%
        function YDiffStructList = stepPeakParameters(obj, options)
            arguments
                obj                    (1,1) Pika.PeakSeries
                options.StepPeakCenter (1,:) {mustBeMember(options.StepPeakCenter, ["Yes", "No"])} = "Yes"
                options.PeakCenterStep (1,1) {mustBeNumeric}                                       = 1
                options.StepPeakWidth  (1,:) {mustBeMember(options.StepPeakWidth, ["Yes", "No"])}  = "Yes"
                options.PeakWidthStep  (1,1) {mustBeNumeric}                                       = 0.01
                options.StepPeakHeight (1,:) {mustBeMember(options.StepPeakHeight, ["Yes", "No"])} = "Yes"
                options.PeakHeightStep (1,1) {mustBeNumeric}                                       = 0.01
            end
            
            import Pika.PeakShape
            
            % Initialize output struct list
            YDiffStructList(1:obj.numPeaks, 1) = struct("peakCenter", [], "peakWidth", [], "peakHeight", []);
            if isempty(obj.peaks)
                return
            end

            % Step each peak and store the output struct in list
            for i = 1:obj.numPeaks
                YDiffStructList(i) = obj.peaks(i).stepPeakParameters("StepPeakCenter", options.StepPeakCenter , ...
                                                                     "PeakCenterStep", options.PeakCenterStep , ...
                                                                     "StepPeakWidth" , options.StepPeakWidth  , ...
                                                                     "PeakWidthStep" , options.PeakWidthStep  , ...
                                                                     "StepPeakHeight", options.StepPeakHeight , ...
                                                                     "PeakHeightStep", options.PeakHeightStep);
            end
        end

        function YDiffStructList = stepSharedPeakParameters(obj, options)
            arguments
                obj                      (1,1) Pika.PeakSeries
                options.StepPeakCenter   (1,:) {mustBeMember(options.StepPeakCenter, ["Yes", "No"])}   = "Yes"
                options.PeakCenterStep   (1,1) {mustBeNumeric}                                         = 1
                options.StepPeakWidth    (1,:) {mustBeMember(options.StepPeakWidth, ["Yes", "No"])}    = "Yes"
                options.PeakWidthStep    (1,1) {mustBeNumeric}                                         = 0.01
                options.StepPeakHeight   (1,:) {mustBeMember(options.StepPeakHeight, ["Yes", "No"])}   = "Yes"
                options.PeakHeightStep   (1,1) {mustBeNumeric}                                         = 0.01
                options.StepPeakSpacing  (1,:) {mustBeMember(options.StepPeakSpacing, ["Yes", "No"])}  = "Yes"
                options.PeakSpacingStep  (1,1) {mustBeNumeric}                                         = 0.01
                options.SharePeakWidth   (1,1) {mustBeMember(options.SharePeakWidth, ["Yes", "No"])}   = "No"
                options.SharePeakSpacing (1,1) {mustBeMember(options.SharePeakSpacing, ["Yes", "No"])} = "No"
            end

            % Initialize output struct list
            YDiffStructList(1:obj.numPeaks, 1) = struct("peakCenter", [], "peakWidth", [], "peakHeight", [], "peakSpacing", []);
            if isempty(obj.peaks)
                return
            end

            % Step each peak and store the output struct in list
            for i = 1:obj.numPeaks
                % Calculate typical parameter steps
                thisYDiffStructList = obj.peaks(i).stepPeakParameters("StepPeakCenter", options.StepPeakCenter , ...
                                                                      "PeakCenterStep", options.PeakCenterStep , ...
                                                                      "StepPeakWidth" , options.StepPeakWidth  , ...
                                                                      "PeakWidthStep" , options.PeakWidthStep  , ...
                                                                      "StepPeakHeight", options.StepPeakHeight , ...
                                                                      "PeakHeightStep", options.PeakHeightStep);
                
                YDiffStructList(i).peakCenter = thisYDiffStructList.peakCenter;
                YDiffStructList(i).peakWidth  = thisYDiffStructList.peakWidth;
                YDiffStructList(i).peakHeight = thisYDiffStructList.peakHeight;

                % Calculate peak spacing step
                if strcmp(options.StepPeakSpacing, "Yes") && strcmp(options.SharePeakSpacing, "Yes")
                    thisPeakCenterStep = (i-1) * options.PeakSpacingStep;
                    peakSpacingYDiffStruct = obj.peaks(i).stepPeakParameters("StepPeakCenter", options.StepPeakSpacing , ...
                                                                             "PeakCenterStep", thisPeakCenterStep      , ...
                                                                             "StepPeakWidth" , "No"                    , ...
                                                                             "StepPeakHeight", "No"                   );
    
                    YDiffStructList(i).peakSpacing = peakSpacingYDiffStruct.peakCenter;
                end
            end
        end

        function obj = stepPeak(obj, parameter, increment, peakIndex)
            arguments
                obj       (1,1) Pika.PeakSeries
                parameter (1,:) {mustBeMember(parameter, ["peakCenter", "peakWidth", "peakHeight"])}
                increment (1,1) {mustBeNumeric}
                peakIndex (1,1) {mustBeInteger, mustBePositive}
            end

            obj.peaks(peakIndex) = obj.peaks(peakIndex).stepPeak(parameter, increment);
        end

        function obj = stepSharedPeaks(obj, parameter, increment)
            arguments
                obj       (1,1) Pika.PeakSeries
                parameter (1,:) {mustBeMember(parameter, ["peakCenter", "peakWidth", "peakHeight"])}
                increment (1,1) {mustBeNumeric}
            end
            
            % Step each peak by the same parameter and increment
            for i = 1:obj.numPeaks
                obj.peaks(i) = obj.peaks(i).stepPeak(parameter, increment);
            end
        end

        %%%---Utility Functions---%%%
        function numPeaks = countPeaks(obj)
            numPeaks = length(obj.peaks);
        end

        function numPts = numPts(obj, options)
            arguments
                obj               (1,1) Pika.Peak
                options.RemoveNaN (1,:) {mustBeMember(options.RemoveNaN, ["Yes", "No"])} = "No"
            end

            numPts = obj.peaks(1).numPts("RemoveNaN", options.RemoveNaN);
        end

        function Y = calculateFit(obj)
            if isempty(obj.peaks)
                Y = [];
                return
            end

            Y = obj.peaks(1).calculateFit();
            for i = 2:length(obj.peaks)
                Y = Y + obj.peaks(i).calculateFit();
            end
        end

        function parameters = getPeakParameters(obj)
            parameters = cell(obj.numPeaks, 5);

            for i = 1:obj.numPeaks
                parameters(i,1) = {i};
                parameters(i, 2:end) = obj.peaks(i).getPeakParameters();
            end
        end

        function spectrumData = getSpectrumData(obj, options)
            arguments
                obj               (1,1) Pika.PeakSeries
                options.RemoveNaN (1,:) {mustBeMember(options.RemoveNaN, ["Yes", "No"])} = "No"
            end

            spectrumData = zeros(obj.peaks(1).numPts("RemoveNaN", options.RemoveNaN), 2 * obj.numPeaks);

            columnIndexPair = [1,2];
            for i = 1:obj.numPeaks
                spectrumData(:, columnIndexPair) = obj.peaks(i).getSpectrumData("RemoveNaN", options.RemoveNaN);
                columnIndexPair = columnIndexPair + 2;
            end
        end
    end
end