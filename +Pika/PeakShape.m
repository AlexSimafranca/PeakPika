classdef PeakShape
    methods (Static)
        function Y = gaussian(X, peakCenter, sigma, options)
            arguments
                X                  (:,1) {mustBeNumeric}
                peakCenter         (1,1) {mustBeNumeric}  = 0
                sigma              (1,1) {mustBeNonnegative} = 0.1
                options.PeakHeight (1,:) {mustBeNumeric}  = 1
            end

            % Calculate Gaussian Y vector
            Y = options.PeakHeight * exp(-(X - peakCenter).^2 / (2 * sigma^2));
        end

        function Y = lorentzian()

        end
    end
end