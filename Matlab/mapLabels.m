function Cx = mapLabels(Cx_truth)
    
    unique_Cx = unique(Cx_truth);
    
    
    label_map_x = containers.Map(unique_Cx, 1:length(unique_Cx));
    
   
    Cx = arrayfun(@(x) label_map_x(x), Cx_truth);
end

