function [cell_id, new_header, cell_data] = prepare_cell_data(input_cell,cell_header)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%cell_id = cell([1,1]);
new_header = cell([1,length(cell_header)-1]);
cell_data = zeros(1,length(cell_header)-1);
jj = 1;
        for i = 1:length(cell_header)
            if strcmp(cell_header{i}, 'filename') == 1 
                cell_id = input_cell.(cell_header{i})(1:end-4);
%             elseif strcmp(cell_header{i},'sweeps_of_interest') == 1
%                 cell_id{2} = num2str(input_cell.(cell_header{i}));
            else              
                cell_data(1,jj) = input_cell.(cell_header{i});
                new_header{jj} = cell_header{i};
                jj = jj + 1;
            end
        end
end

