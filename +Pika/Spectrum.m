classdef Spectrum
    properties (SetAccess = private)
        X      (:,1) double
        Y      (:,1) double
        ID     (1,1) int8   
        name   (1,:) string 
        XUnits (1,:) string 
        YUnits (1,:) string 
        colors (:,1) cell
    end

    methods (Static)
        function obj = Spectrum(X, Y, name, ID, options)
            arguments
                X              (:,1) {mustBeNumeric} = []
                Y              (:,1) {mustBeNumeric} = []
                name           (1,:) string          = ""
                ID             (1,1) int8            = 0
                options.XUnits (1,:) string          = "nm"
                options.YUnits (1,:) string          = "a.u."
            end

            obj.X = X;
            obj.Y = Y;
            obj.name = name;
            obj.ID = ID;
            obj.XUnits = options.XUnits;
            obj.YUnits = options.YUnits;
            
            % If X provided but no Y, ensure that Y matches X in length
            if isempty(Y) && ~isempty(X)
                obj.Y = NaN(length(X), 1);
            end
            
            % Set default plot colors
            obj.colors = {[247, 150,  70] ./ 255 ; ...
                          [242, 207,   0] ./ 255 ; ...
                          [121, 192,  40] ./ 255 ; ...
                          [  0, 176,  80] ./ 255 ; ... 
                          [ 79, 129, 189] ./ 255 ; ...
                          [128, 100, 162] ./ 255 ; ...
                          [ 50,  50,  50] ./ 255};
        end

        %%%---File Reading Functions---%%%
        function obj = read(filePath, options)
            arguments
                filePath                (1,:) {mustBeText}
                options.NamePrelimiter  (1,:) {mustBeText}    = ""
                options.NamePostlimiter (1,:) {mustBeText}    = ""
                options.HeaderLines     (1,1) {mustBeInteger} = 0
            end

            import Pika.Spectrum

            % Reformat data name
            dataName = Spectrum.formatName(filePath, options.NamePrelimiter, options.NamePostlimiter);

            % Read data from file
            data = readmatrix(filePath, 'NumHeaderLines', options.HeaderLines);

            % Partition data into spectra
            dataStruct = Spectrum.partitionData(data);
            numSpectra = length(fields(dataStruct));

            obj(numSpectra, 1) = Spectrum();
            for i = 1:numSpectra
                thisData = dataStruct.("spectrum" + string(i));

                % Create new spectrum object
                newObj = Spectrum();
                newObj.X = thisData(:,1);
                newObj.Y = thisData(:,2);
                newObj.name = dataName + "\_" + string(i);

                % Append new spectrum to output array
                obj(i) = newObj;
            end

            % If only one spectrum, remove suffix from name
            if numSpectra == 1
                obj(1).name = dataName;
            end
        end

        function name = formatName(filePath, prelimiter, postlimiter)
            arguments
                filePath    (1,:) {mustBeText}
                prelimiter  (1,:) {mustBeText} = ""
                postlimiter (1,:) {mustBeText} = ""
            end
            
            % Extract name from file path
            [~, rawName, ~] = fileparts(filePath);

            % Trim name to portion between limiters
            if ~strcmp(prelimiter, "") && contains(rawName, prelimiter)
                rawName = extractAfter(rawName, prelimiter);
            end

            if ~strcmp(postlimiter, "") && contains(rawName, postlimiter)
                rawName = extractBefore(rawName, postlimiter);
            end

            % Escape underscores to prevent formatting issues
            name = string(replace(rawName, "_", "\_"));
        end

        function dataStruct = partitionData(dataMatrix)
            arguments
                dataMatrix (:,:) {mustBeNumeric}
            end

            % Return early if data already single spectrum
            numCols = size(dataMatrix, 2);
            if numCols == 2
                dataStruct = struct("spectrum1", dataMatrix);
                return
            end

            % Check that data matrix is partitionable
            assert(mod(numCols, 2) == 0, ...
                   "[ Spectrum > partitionData ] Data partition error: Data " + ...
                   "matrix cannot be partitioned into spectra because it "    + ...
                   "contains an odd number of columns.")

            % Separate data in pairs of columns
            dataStruct = struct();
            for i = 1:numCols / 2
                dataStruct.("spectrum" + string(i)) = dataMatrix(:, 2 * (i-1) + 1 : 2*i);
            end
        end
    
        %%%---Utility Functions---%%%
        function [rangeMask1, rangeMask2] = findCommonRange(X1, X2)
            arguments
                X1 (1,:) {mustBeNumeric}
                X2 (1,:) {mustBeNumeric}
            end
            
            import Pika.Spectrum

            numPts1 = length(X1);
            numPts2 = length(X2);
            rangeMask1 = false(numPts1, 1);
            rangeMask2 = false(numPts2, 1);

            point1 = 1;
            point2 = 1;
            while point1 <= numPts1 && point2 <= numPts2
                XValue1 = X1(point1);
                XValue2 = X2(point2);

                % Check for point match
                if abs(XValue1 - XValue2) < 1E4 * min(eps(XValue1), eps(XValue2))
                    % Indicate match in range masks
                    rangeMask1(point1) = true;
                    rangeMask2(point2) = true;
                    
                    % Advance both pointers
                    point1 = point1 + 1;
                    point2 = point2 + 1;
                    continue
                end

                % Advance pointer for position with lower value
                if XValue1 < XValue2
                    point1 = point1 + 1;
                else
                    point2 = point2 + 1;
                end
            end
        end

        function index = findXValue(X, value)
            arguments
                X                 (:,1) {mustBeNumeric}
                value             (1,1) {mustBeNumeric}
            end

            diffVect = abs(X - value);
            index = find(diffVect == min(diffVect), 1);
        end

        function axesHandle = setAxes(axesHandle, options)
            arguments
                axesHandle        (1,1) matlab.graphics.axis.Axes                        = matlab.graphics.axis.Axes()
                options.NewFigure (1,:) {mustBeMember(options.NewFigure, ["Yes", "No"])} = "Yes"
            end

            % If no axes handle passed, reference current or create new axes
            if isempty(axesHandle.Parent)
                if strcmp(options.NewFigure, "Yes")
                    figure()
                end
                
                axesHandle = gca;
            end
        end
    
        function color = getColor(colorIndex)
            import Pika.Spectrum

            dummySpectrum = Spectrum();
            numColors = length(dummySpectrum.colors);
            colorIndex = mod(colorIndex, numColors) + 1;

            color = dummySpectrum.colors{colorIndex};
        end
    end

    methods 
        %%%---Property Access Functions---%%%
        % Spectrum identifiers
        function obj = setID(obj, ID)
            arguments
                obj (1,1) Pika.Spectrum
                ID  (1,1) {mustBeInteger}
            end

            obj.ID = uint8(ID);
        end

        function obj = setName(obj, name)
            arguments
                obj  (1,1) Pika.Spectrum
                name (1,:) {mustBeText}
            end

            obj.name = name;
        end

        % X and Y data
        function obj = setX(obj, X)
            arguments
                obj (1,1) Pika.Spectrum
                X   (:,1) {mustBeNumeric}
            end

            assert(length(X) == length(obj.Y), ...
                   "[ Spectrum > setXY ] Data setting error: X and Y vectors must "     + ...
                   "have the same length. To set both X and Y vectors simultaneously, " + ...
                   "use the setXY() method.")

            obj.X = X;
        end

        function obj = setY(obj, Y)
            arguments
                obj (1,1) Pika.Spectrum
                Y   (:,1) {mustBeNumeric}
            end

            % If empty Y vector passed, clear content of Y vector
            if isempty(Y)
                obj.Y = NaN(length(obj.X), 1);
                return
            end

            % Otherwise, replace Y data with passed Y vector
            assert(length(obj.X) == length(Y), ...
                   "[ Spectrum > setXY ] Data setting error: X and Y vectors must "     + ...
                   "have the same length. To set both X and Y vectors simultaneously, " + ...
                   "use the setXY() method.")

            obj.Y = Y;
        end

        function obj = updateY(obj, Y, selectionMask)
            arguments
                obj           (1,1) Pika.Spectrum
                Y             (:,1) {mustBeNumeric}
                selectionMask (:,1) logical         = true(0,1)
            end

            if ~isempty(selectionMask)
                % Check that selection mask is equal in length to Y data
                assert(length(obj.Y) == length(selectionMask), ...
                       "[ Spectrum > updateY ] Y data update error: Selection mask " + ...
                       "must be the same length as the stored Y data vector.")

                % Check that input Y vector has enough points to map with selection vector
                numSelectedPoints = sum(selectionMask);
                assert(length(Y) == numSelectedPoints, ...
                       "[ Spectrum > updateY ] Y data update error: Number of new data " + ...
                       "points in input Y vector must equal the number of data points "  + ...
                       "selected by the selection mask (" + string(numSelectedPoints) + ")")
                
                obj.Y(selectionMask) = Y;
                return
            end

            obj = obj.setY(Y);
        end

        function obj = setXY(obj, X, Y)
            arguments
                obj (1,1) Pika.Spectrum
                X   (:,1) {mustBeNumeric}
                Y   (:,1) {mustBeNumeric}
            end

            assert(length(X) == length(Y), ...
                   "[ Spectrum > setXY ] Data setting error: X and Y vectors must " + ...
                   "have the same length.")

            obj.X = X;
            obj.Y = Y;
        end

        %%%---Data Manipulation Functions---%%%
        function obj = normalize(obj, maxRange, options)
            arguments 
                obj                  (1,1) Pika.Spectrum
                maxRange             (1,2) {mustBeNumeric}                                     = [NaN, NaN]
                options.MinRange     (1,2) {mustBeNumeric}                                     = [NaN, NaN]
                options.SetMinToZero (1,:) {mustBeMember(options.SetMinToZero, ["Yes", "No"])} = "No"
                options.NameTag      (1,:) {mustBeText}                                        = "\_norm"
            end

            % If placeholder ranges, replace with full spectrum range
            if any(isnan(maxRange), "all")
                maxRange = [min(obj.X), max(obj.X)];
            end

            if any(isnan(options.MinRange), "all")
                options.MinRange = [min(obj.X), max(obj.X)];
            end

            % Calculate max and min values within specified ranges
            maxMask = obj.X >= maxRange(1) & obj.X <= maxRange(2);
            minMask = obj.X >= options.MinRange(1) & obj.X <= options.MinRange(2);
            YMax = max(obj.Y(maxMask));
            YMin = min(obj.Y(minMask));

            % Zero min value, then normalize to max value
            if strcmp(options.SetMinToZero, "Yes")
                obj.Y = obj.Y - YMin;
                obj.Y = obj.Y / (YMax - YMin);
                return
            end

            obj.Y = obj.Y / YMax;
            obj.name = obj.name + options.NameTag;
        end

        function obj = trimData(obj, dataRange, options)
            arguments
                obj               (1,1) Pika.Spectrum
                dataRange         (1,2) {mustBeNumeric}
                options.RangeType (1,:) {mustBeMember(options.RangeType, ["Value", "Index"])} = "Value"
            end
            
            % If indices provided, convert to values
            if strcmp(options.RangeType, "Index")
                options.RangeType = [obj.X(dataRange(1)), obj.X(dataRange(2))];
            end
            
            % Check validity of specified data range
            assert(dataRange(2) ~= dataRange(1), ...
                   "[ Spectrum > trimData ] Data trimming error: Specified data range " + ...
                   "limits are identical: [" + string(dataRange(1)) + ", "              + ...
                   string(dataRange(2)) + "]")

            assert(dataRange(2) > dataRange(1), ...
                   "[ Spectrum > trimData ] Data trimming error: Specified data range "   + ...
                   "limits [" + string(dataRange(1)) + ", " + string(dataRange(2)) + "] " + ...
                   "are not in ascending order.")

            % Set data outside of range to NaN
            isInRange = obj.X >= dataRange(1) & obj.X <= dataRange(2);
            obj.X(~isInRange) = NaN;
            obj.Y(~isInRange) = NaN;

            % Append suffix to name
            obj.name = obj.name + "\_trim";
        end

        function obj = maskData(obj, maskRange)
            arguments
                obj       (1,1) Pika.Spectrum
                maskRange (1,2) {mustBeNumeric}
            end

            isInRange = obj.X > maskRange(1) & obj.X < maskRange(2);
            obj.Y(isInRange) = NaN;
        end

        % X-axis conversions
        function obj = convert2eV(obj)
            import Pika.Units

            if strcmp(obj.XUnits, "eV")
                disp("[ Spectrum > convert2eV ] X-axis conversion error: X-axis already in units of eV.")
                return
            end

            obj.X = Units.nm2eV(obj.X);
            obj.XUnits = "eV";
        end

        function obj = convert2nm(obj)
            import Pika.Units

            if strcmp(obj.XUnits, "nm")
                disp("[ Spectrum > convert2eV ] X-axis conversion error: X-axis already in units of nm.")
                return
            end

            obj.X = Units.eV2nm(obj.X);
            obj.XUnits = "nm";
        end

        %%%---Plotting Functions---%%%
        function plotHandle = plot(obj, axesHandle, options)
            arguments
                obj               (1,1) Pika.Spectrum
                axesHandle        (1,1) matlab.graphics.axis.Axes                        = matlab.graphics.axis.Axes()
                options.NewFigure (1,:) {mustBeMember(options.NewFigure, ["Yes", "No"])} = "Yes"
            end

            import Pika.Spectrum

            % Plot data
            axesHandle = Spectrum.setAxes(axesHandle, "NewFigure", options.NewFigure);
            plotHandle = plot(axesHandle, obj.X, obj.Y, "DisplayName", obj.name , "LineWidth", 3);
            hold(axesHandle, "On")

            % Apply plot settings
            switch obj.XUnits
                case "nm"
                    xlabel(axesHandle, "Wavelength (nm)", "FontWeight", "Bold")
                case "eV"
                    xlabel(axesHandle, "Energy (eV)", "FontWeight", "Bold")
            end
            
            ylabel(axesHandle, "Absorbance (" + obj.YUnits + ")", "FontWeight", "Bold")

            set(axesHandle, "Box", "on", "LineWidth", 1.5, "FontSize", 20)
            legend(axesHandle, "Show")
        end
    
        %%%---Utility Functions---%%%
        function bool = isEmpty(obj)
            bool = isempty(obj.X) && isempty(obj.Y);
        end
    
        function numPts = numPts(obj, options)
            arguments
                obj               (1,1) Pika.Spectrum
                options.RemoveNaN (1,:) {mustBeMember(options.RemoveNaN, ["Yes", "No"])} = "No"
            end            
            
            if strcmp(options.RemoveNaN, "Yes")
                isNanData = isnan(obj.X) & isnan(obj.Y);
                numPts = sum(~isNanData);
            else
                numPts = length(obj.X);
            end
        end

        function spectrumData = getSpectrumData(obj, options)
            arguments
                obj                (1,1) Pika.Spectrum
                options.ReturnType (1,:) {mustBeMember(options.ReturnType, ["Matrix", "Table"])} = "Matrix"
                options.RemoveNaN  (1,:) {mustBeMember(options.RemoveNaN, ["Yes", "No"])} = "No"
            end

            spectrumData = [obj.X, obj.Y];

            if strcmp(options.RemoveNaN, "Yes")
                isNanData = isnan(obj.X) & isnan(obj.Y);
                spectrumData = spectrumData(~isNanData, :);
            end

            % Convert to table, if requested
            if strcmp(options.ReturnType, "Table")
                dataTable = table();
                dataTable(:,1) = table(spectrumData(:,1));
                dataTable(:,2) = table(spectrumData(:,2));
                dataTable.Properties.VariableNames = ["X", "Y"];
    
                spectrumData = dataTable;
            end
        end
    end
end