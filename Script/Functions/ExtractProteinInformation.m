
%% Extract protein information

function infoTable=ExtractProteinInformation(strs)


info = struct('ID',{},'locus_tag', {}, 'protein', {}, 'protein_id', {}, 'location', {},'gene',{});

for i = 1:length(strs)
    str = strs{i};
    
    spacePos = find(str == ' ', 1);
if ~isempty(spacePos)
    extractedPart = str(1:spacePos-1);
else
    extractedPart = 'nan';
end

    locus_tag = regexp(str, '\[locus_tag=(.*?)\]', 'tokens', 'once');
    if ~isempty(locus_tag)
        locus_tag = locus_tag{1};
    else
        locus_tag = 'nan';
    end
    
 
    protein = regexp(str, '\[protein=(.*?)\]', 'tokens', 'once');
    if ~isempty(protein)
        protein = protein{1};
    else
        protein = 'nan';
    end
  
    protein_id = regexp(str, '\[protein_id=(.*?)\]', 'tokens', 'once');
    if ~isempty(protein_id)
        protein_id = protein_id{1};
    else
        protein_id = 'nan';
    end

    location = regexp(str, '\[location=(.*?)\]', 'tokens', 'once');
    if ~isempty(location)
        location = location{1};
    else
        location = 'nan';
    end    
    

    gene = regexp(str, '\[gene=(.*?)\]', 'tokens', 'once');
    if ~isempty(gene)
        gene = gene{1};
    else
        gene ='nan';
    end


    info(i).ID = extractedPart;
    info(i).locus_tag = locus_tag;
    info(i).protein = protein;
    info(i).protein_id = protein_id;
    info(i).location = location;
    info(i).gene = gene;
    
    infoTable = struct2table(info);
end


disp(info);
end
