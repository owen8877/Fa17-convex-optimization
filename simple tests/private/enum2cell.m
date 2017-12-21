function[out] = enum2cell(e)
% 'Unpack' a Java enumeration into a Matlab cell array
% Example : none (Oh, the FEX code metrics..) 
    c = {};
    while e.hasMoreElements
        c{end+1} = e.nextElement;                %#ok 
    end 
    out = c(:);
end    
