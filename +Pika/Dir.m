classdef (Abstract) Dir
    methods (Static)
        function filenames = read(path, fileExt)
            arguments
                path    (1,:) {mustBeText}
                fileExt (1,:) {mustBeText} = [".txt", ".csv", ".dat", ".xls", ".xlsx", ".xlsm"]
            end

            % Query directory for files with given extension
            filenames = strings(100,1);
            i = 1;
            for ext = fileExt
                % Retrieve files with extension in directory
                pathQuery = path + "/*" + ext;
                dirContent = dir(pathQuery);

                if ~isempty(dirContent)
                    % Store filenames in a cell array
                    newFilenames = cellfun(@(x) string(x), {dirContent.name});
                    numNewNames = length(newFilenames);
    
                    % Append new filenames to output array
                    filenames(i : i + numNewNames-1) = newFilenames;
                    i = i + numNewNames;
                end
            end

            % Remove empty preallocated values
            isEmptyStr = strcmp(filenames, "");
            filenames(isEmptyStr) = [];
        end

        function [sortedlist, dataname] = namesort(namelist, prelim, postlim)
            digits = ['1','2','3','4','5','6','7','8','9','0'];
            
            % Extract dataset names
            n_items = length(namelist);
            dataname = cell(1, n_items);
            
            if ~isempty(prelim) && ~isempty(postlim)
                for i = 1:n_items
                    dataname{i} = extractBefore(namelist{i}, postlim);
                    dataname{i} = extractAfter(dataname{i}, prelim);
                end
            elseif ~isempty(prelim)
                for i = 1:n_items
                    dataname{i} = extractAfter(dataname{i}, prelim);
                end
            elseif ~isempty(postlim)
                for i = 1:n_items
                    dataname{i} = extractBefore(namelist{i}, postlim);
                end
            else
                dataname = namelist;
            end

            % Insert a zero where only single numerical digit
            for i = 1:n_items
                % Count number of digits in sample number
                isdigit = ismember(dataname{i}, digits);  % check for numerical characters
                digitcount = sum(isdigit);
                
                if digitcount == 1
                    % Insert a zero in front of single-digit number
                    index = 1:length(isdigit);
                    n_digit = index(isdigit);
            
                    dataname{i} = [dataname{i}(1:n_digit-1), '0', dataname{i}(n_digit:end)];
                end
            end
            
            % Store namelist and dataname together for sorting
            list_temp = [dataname', namelist'];
            
            % Sort temporary list
            list_temp = sortrows(list_temp);
            
            % Pass sorted lists to output variables
            dataname = cell(n_items, 1);
            sortedlist = cell(n_items, 1);
            
            for i = 1:n_items
                dataname{i} = list_temp{i};
                sortedlist{i} = list_temp{n_items + i};
            end
        end

        function filename = getFilenameFromPath(path)
            arguments
                path (1,:) {mustBeText}
            end

            % Split path into directory and file components
            [~, file, ext] = fileparts(path);
            assert(~strcmp(file, "") && ~strcmp(ext, ""), ...
                   "[ Dir > getFilenameFromPath ] Filename retrieval error: " + ...
                   "Path does not contain a filename.")

            % Reconstruct filename
            filename = file + ext;
        end
    end
end