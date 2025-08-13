classdef Peak
    properties
        peakSpectrum (1,1) Pika.Spectrum
        peakCenter   (1,1) {mustBeNumeric}
        sigma        (1,1) {mustBeNumeric, mustBeNonnegative}
        peakHeight   (1,:) {mustBeNumeric}  
        peakType     (1,:) string
    end

    methods (Static)
        function obj = Peak(X, peakCenter, sigma, options)
            arguments
                X                  (:,1) {mustBeNumeric}                                              = []
                peakCenter         (1,1) {mustBeNumeric}                                              = 0
                sigma              (1,1) {mustBeNumeric, mustBeNonnegative}                           = 1
                options.PeakType   (1,:) {mustBeMember(options.PeakType, ["Gaussian", "Lorentzian"])} = "Gaussian"
                options.PeakHeight (1,:) {mustBeNumeric}                                              = 1
            end

            import Pika.Spectrum
            import Pika.PeakShape

            % Store peak parameters on object
            obj.peakCenter = peakCenter;
            obj.sigma = sigma;
            obj.peakHeight = options.PeakHeight;
            obj.peakType = options.PeakType;

            % Generate peak spectrum depending on peak type
            if isempty(X)
                obj.peakSpectrum = Spectrum(X);
                return
            end

            obj.peakType = options.PeakType;
            switch options.PeakType
                case "Gaussian"
                    obj.peakSpectrum = Spectrum(X, [], "Gaussian");

                    % Calculate gaussian Y data
                    Y = PeakShape.gaussian(X, peakCenter, sigma, "PeakHeight", options.PeakHeight);
                    obj.peakSpectrum = obj.peakSpectrum.setY(Y);
                case "Lorentzian"
                    obj.peakSpectrum = Spectrum(X);
            end
        end
    end

    methods
        %%%---Data Manipulation Functions---%%%
        function obj = recalculatePeak(obj)
            arguments
                obj (1,1) Pika.Peak
            end

            import Pika.Spectrum
            import Pika.PeakShape

            % Recalculate peak Y data
            X = obj.peakSpectrum.X;
            switch obj.peakType
                case "Gaussian"
                    YNew = PeakShape.gaussian(X, obj.peakCenter, obj.sigma, "PeakHeight", obj.peakHeight);
                    obj.peakSpectrum = obj.peakSpectrum.setY(YNew);
                case "Lorentzian"
                    obj.peakSpectrum = Spectrum(X);
            end
        end

        % X-axis conversions
        function obj = convert2eV(obj)
            obj.peakSpectrum = obj.peakSpectrum.convert2eV();
        end

        function obj = convert2nm(obj)
            obj.peakSpectrum = obj.peakSpectrum.convert2nm();
        end

        %%%---Plotting Functions---%%%
        function plotHandle = plot(obj, axesHandle, options)
            arguments
                obj               (1,1) Pika.Peak
                axesHandle        (1,1) matlab.graphics.axis.Axes                        = matlab.graphics.axis.Axes()
                options.NewFigure (1,:) {mustBeMember(options.NewFigure, ["Yes", "No"])} = "Yes"
            end

            import Pika.Spectrum

            axesHandle = Spectrum.setAxes(axesHandle, "NewFigure", options.NewFigure);
            plotHandle = obj.peakSpectrum.plot(axesHandle);
        end

        %%%---Peak Fitting Functions---%%%
        function YStruct = stepPeakParameters(obj, options)
            arguments
                obj                    (1,1) Pika.Peak
                options.StepPeakCenter (1,:) {mustBeMember(options.StepPeakCenter, ["Yes", "No"])} = "Yes"
                options.PeakCenterStep (1,:) {mustBeNumeric}                                       = 1
                options.StepPeakWidth  (1,:) {mustBeMember(options.StepPeakWidth, ["Yes", "No"])}  = "Yes"
                options.PeakWidthStep  (1,:) {mustBeNumeric}                                       = 0.01
                options.StepPeakHeight (1,:) {mustBeMember(options.StepPeakHeight, ["Yes", "No"])} = "Yes"
                options.PeakHeightStep (1,:) {mustBeNumeric}                                       = 0.01
            end
            
            import Pika.PeakShape
            
            X = obj.peakSpectrum.X;
            YStruct = struct("peakCenter", [], "peakWidth", [], "peakHeight", []);

            % Select peak shape
            switch obj.peakType
                case "Gaussian"
                    peakFunction = @PeakShape.gaussian;
                case "Lorentzian"
                    peakFunction = @PeakShape.lorentzian;
            end

            % Step peak parameters and store Y data difference
            if strcmp(options.StepPeakCenter, "Yes")
                % If step size passed, set symmetric steps
                if length(options.PeakCenterStep) ~= 2
                    options.PeakCenterStep = [-options.PeakCenterStep(1), options.PeakCenterStep(1)];
                end

                % Apply steps
                YPlus = peakFunction(X, obj.peakCenter + options.PeakCenterStep(2), obj.sigma, ...
                                     "PeakHeight", obj.peakHeight);
                YMinus = peakFunction(X, obj.peakCenter + options.PeakCenterStep(1), obj.sigma, ...
                                      "PeakHeight", obj.peakHeight); 
                YStruct.("peakCenter") = [YPlus - obj.peakSpectrum.Y, YMinus - obj.peakSpectrum.Y];
            end
            
            if strcmp(options.StepPeakWidth, "Yes")
                % If step size passed, set symmetric steps
                if length(options.PeakWidthStep) ~= 2
                    options.PeakWidthStep = [-options.PeakWidthStep(1), options.PeakWidthStep(1)];
                end

                % Apply steps
                YPlus = peakFunction(X, obj.peakCenter, obj.sigma + options.PeakWidthStep(2), ...
                                     "PeakHeight", obj.peakHeight);
                YMinus = peakFunction(X, obj.peakCenter, max(obj.sigma + options.PeakWidthStep(1), abs(options.PeakWidthStep(1))), ...
                                      "PeakHeight", obj.peakHeight); 
                YStruct.("peakWidth") = [YPlus - obj.peakSpectrum.Y, YMinus - obj.peakSpectrum.Y];
            end

            if strcmp(options.StepPeakHeight, "Yes")
                % If step size passed, set symmetric steps
                if length(options.PeakHeightStep) ~= 2
                    options.PeakHeightStep = [-options.PeakHeightStep(1), options.PeakHeightStep(1)];
                end

                % Apply steps
                YPlus = peakFunction(X, obj.peakCenter, obj.sigma, ...
                                     "PeakHeight", obj.peakHeight + options.PeakHeightStep(2));
                YMinus = peakFunction(X, obj.peakCenter, obj.sigma, ...
                                      "PeakHeight", obj.peakHeight + options.PeakHeightStep(1)); 
                YStruct.("peakHeight") = [YPlus - obj.peakSpectrum.Y, YMinus - obj.peakSpectrum.Y];
            end
        end

        function obj = stepPeak(obj, parameter, increment)
            arguments
                obj       (1,1) Pika.Peak
                parameter (1,:) {mustBeMember(parameter, ["peakCenter", "peakWidth", "peakHeight"])}
                increment (1,1) {mustBeNumeric}
            end
         
            switch parameter
                case "peakCenter"
                    obj.peakCenter = obj.peakCenter + increment;
                case "peakWidth"
                    obj.sigma = obj.sigma + increment;
                case "peakHeight"
                    obj.peakHeight = obj.peakHeight + increment;
            end

            obj = obj.recalculatePeak();
        end

        %%%---Utility Functions---%%%
        function numPeaks = countPeaks(obj)
            numPeaks = length(obj.peakSpectrum);
        end

        function numPts = numPts(obj, options)
            arguments
                obj               (1,1) Pika.Peak
                options.RemoveNaN (1,:) {mustBeMember(options.RemoveNaN, ["Yes", "No"])} = "No"
            end

            numPts = obj.peakSpectrum.numPts("RemoveNaN", options.RemoveNaN);
        end

        function Y = calculateFit(obj)
            Y = obj.peakSpectrum.Y;
        end

        function parameters = getPeakParameters(obj)
            parameters = {obj.peakType, obj.peakCenter, obj.sigma, obj.peakHeight};
        end

        function spectrumData = getSpectrumData(obj, options)
            arguments
                obj               (1,1) Pika.Peak
                options.RemoveNaN (1,:) {mustBeMember(options.RemoveNaN, ["Yes", "No"])} = "No"
            end

            spectrumData = obj.peakSpectrum.getSpectrumData("RemoveNaN", options.RemoveNaN);
        end
    end
end