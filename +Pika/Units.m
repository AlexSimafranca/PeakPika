classdef (Abstract) Units
    properties (Constant)
        % Physical constants
        h = 6.62607015E-34;   % m^2*kg/s
        c = 2.99792458E8;   % m/s

        % Conversion constants
        joulePerEV = 1.602176634E-19;   % J/eV
    end

    methods (Static)
        %%%---Unit Conversions---%%%
        %---Energy Conversions---%
        % eV-Joules
        function eV = joule2eV(joule)
            arguments
                joule (:,:) {mustBeNumeric}
            end

            import Pika.Units

            eV = Units.joulePerEV * joule;
        end

        function joule = eV2joule(eV)
            arguments
                eV (:,:) {mustBeNumeric}
            end

            import Pika.Units

            joule = eV / Units.joulePerEV;
        end

        %---Photon Conversions---%
        % nm-eV
        function energy = nm2eV(wavelength)
            arguments 
                wavelength (:,:) {mustBeNumeric}
            end

            import Pika.Units

            energy = (Units.c * Units.h * 1E9 / Units.joulePerEV) ./ wavelength;
        end

        function wavelength = eV2nm(energy)
            arguments 
                energy (:,:) {mustBeNumeric}
            end

            import Pika.Units

            wavelength = (Units.c * Units.h * 1E9 / Units.joulePerEV) ./ energy;
        end
        
        % invcm-eV
        function energy = invcm2eV(frequency)
            arguments
                frequency (:,:) {mustBeNumeric}
            end
                
            import Pika.Units

            energy = (Units.h * Units.c * 100 / Units.joulePerEV) * frequency;
        end

        function frequency = eV2invcm(energy)
            arguments
                energy (:,:) {mustBeNumeric}
            end
                
            import Pika.Units

            frequency = energy / (Units.h * Units.c * 100 / Units.joulePerEV);
        end
    end
end