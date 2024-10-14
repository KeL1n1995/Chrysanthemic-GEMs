
%% Extract protein information
% 输入字符串
function infoTable=ExtractProteinInformation(strs)

% strs 为基因组中的注释信息
% 初始化存储提取信息的结构体数组
info = struct('ID',{},'locus_tag', {}, 'protein', {}, 'protein_id', {}, 'location', {},'gene',{});

% 提取信息
for i = 1:length(strs)
    str = strs{i};
    
    spacePos = find(str == ' ', 1);
if ~isempty(spacePos)
    extractedPart = str(1:spacePos-1);
else
    extractedPart = 'nan';
end


    % 提取 locus_tag
    locus_tag = regexp(str, '\[locus_tag=(.*?)\]', 'tokens', 'once');
    if ~isempty(locus_tag)
        locus_tag = locus_tag{1};
    else
        locus_tag = 'nan';
    end
    
    % 提取 protein
    protein = regexp(str, '\[protein=(.*?)\]', 'tokens', 'once');
    if ~isempty(protein)
        protein = protein{1};
    else
        protein = 'nan';
    end
    
    % 提取 protein_id
    protein_id = regexp(str, '\[protein_id=(.*?)\]', 'tokens', 'once');
    if ~isempty(protein_id)
        protein_id = protein_id{1};
    else
        protein_id = 'nan';
    end
    
    % 提取 location
    location = regexp(str, '\[location=(.*?)\]', 'tokens', 'once');
    if ~isempty(location)
        location = location{1};
    else
        location = 'nan';
    end    
    
    % 提取 gene id
    gene = regexp(str, '\[gene=(.*?)\]', 'tokens', 'once');
    if ~isempty(gene)
        gene = gene{1};
    else
        gene ='nan';
    end
    
    % 存储到结构体数组中
    info(i).ID = extractedPart;
    info(i).locus_tag = locus_tag;
    info(i).protein = protein;
    info(i).protein_id = protein_id;
    info(i).location = location;
    info(i).gene = gene;
    
    infoTable = struct2table(info);
end

% 显示提取的信息
disp(info);
end
