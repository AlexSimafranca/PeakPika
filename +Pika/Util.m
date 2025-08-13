classdef (Abstract) Util
    methods (Static)
        function [newPlotHandleList, newAxesHandle] = transferPlotHandles(plotHandleList, newAxesHandle)
            arguments
                plotHandleList (1,:) matlab.graphics.Graphics
                newAxesHandle  (1,1) matlab.graphics.axis.Axes = matlab.graphics.axis.Axes()
            end

            % If axes handle has no parent figure, generate new figure
            if isempty(newAxesHandle.Parent)
                figure()
                newAxesHandle = gca;
            end

            % Duplicate plots referenced by handles list in new axes
            numPlotHandles = length(plotHandleList);
            newPlotHandleList = gobjects(numPlotHandles, 1);

            for i = 1:numPlotHandles
                newPlotHandleList(i) = copyobj(plotHandleList(i), newAxesHandle);
            end
        end
    end
end