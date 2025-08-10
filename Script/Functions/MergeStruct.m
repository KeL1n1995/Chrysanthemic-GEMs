function mergedResults=MergeStruct(results)

fields = fieldnames(results);
mergedResults = {};
for i = 1:numel(fields)
    fieldName = fields{i};
    mergedResults = [mergedResults ; results.(fieldName)];
end
end
