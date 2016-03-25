function result = contains(mat, val)

result = ismember(mat, val);

for i = 1:ndims(mat)
    result = any(result);
end
  