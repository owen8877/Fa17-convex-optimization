classdef myHashtable

    properties
        h
    end

    methods

        function[out] = myHashtable(varargin)
            % Construct a hash table, optionally with an initial set of keys and values
            % each supplied in a cell array
            
            assert(ismember(nargin,[0 2]),'??? Invalid naumber of arguments.') 
            out.h = java.util.Hashtable;
         
            if nargin > 0
               keys = varargin{1};
               vals = varargin{2};
                  
               assert(iscell(keys),'??? ''keys'' must be supplied in a cell array.')
               assert(iscell(vals),'??? ''values'' must be supplied in a cell array.')
                  
               n = length(keys);
               assert(length(vals) == n,'??? ''keys''and ''values'' must have the same length.')
              
               for i = 1:n
                   out.h.put(keys{i},vals{i});
               end
               
            end
        end

        function clear(obj)
            % Clear hash table
            obj.h.clear;
        end

        function[out] = fieldnames(obj)          %#ok method overloaded to hide h
            out = {};
        end

        function[out] = iskey(obj,keys)
            % Test if hash table maps any values to specified objects
            if iscell(keys)
                fun = @(x) obj.h.containsKey(x);
                out = get(obj,keys,fun);
                out = logical(cell2mat(out));
            else
                out = obj.h.containsKey(keys);
            end
        end

        function[out] = isvalue(obj,values)
            % Test if hash table maps specified objects to any keys
           if iscell(values)
                fun = @(x) obj.h.containsValue(x);
                out = get(obj,values,fun);
                out = logical(cell2mat(out));
            else
                out = obj.h.containsValue(values);
            end
        end

        function[out] = get(obj,keys)
            % Return mapped-to values of specified keys
            fun = @(x) obj.h.get(x);
            out = get(obj,keys,fun);
        end

        function put(obj,keys,values)
            % Map specified key to specified value
            if iscell(keys)
               for i = 1:length(keys)
                   obj.h.put(keys{i},values{i});
               end
            else
               obj.h.put(keys,values);
            end
        end

        function[out] = isempty(obj)
            % Test if hash table is empty
            out = obj.h.isEmpty;
        end

        function[out] = keys(obj)
            % Return hash table's keys
            out = enum2cell(obj.h.keys);
        end

        function rehash(obj)
            % Reorganize hash table, to access its entries more efficiently
            obj.h.rehash;
        end

        function remove(obj,keys)
            % Remove a key (and corresponding value) from hash table
            if iscell(keys)
               for i = 1:length(keys)
                   obj.h.remove(keys{i});
               end
            else
               obj.h.remove(keys);
            end
        end

        function[out] = size(obj)
            % Return the number of keys in hash table
            out = obj.h.size;
        end

        function[out] = values(obj)
            % Return hash table's values
            out = enum2cell(obj.h.elements);
        end

    end

end