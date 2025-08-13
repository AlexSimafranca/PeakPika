classdef (Abstract) Stats
    methods (Static)
        function RMSD = calculateRMSD(data1, data2)
            arguments
                data1 (:,:) {mustBeNumeric}
                data2 (:,:) {mustBeNumeric}
            end

            % Check size of provided data to compare
            assert(all(size(data1) == size(data2)), ...
                   "[ Stats > variance ] Variance calculation error: Data to compare are not the same size.")
                        
            % Calculate difference at each data point
            squareDiff = (data2 - data1).^2;

            % Remove NaN data points from calculation
            squareDiff(isnan(squareDiff)) = [];

            % Calculate root-mean-square difference
            RMSD = sqrt(mean(mean(squareDiff)));
        end
    end
end