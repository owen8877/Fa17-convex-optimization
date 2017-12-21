function[out] = get(obj,x,fun)                   %#ok
if iscell(x)
   out = cellfun(fun,x,'UniformOutput',false);
   out = out(:);
else
   out = fun(x);
end