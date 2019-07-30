import pandas as pd

def from_GenomicCoord_list_to_DataFrame(gen_coord_list):
    chr_list = []
    start_list = []
    end_list = []
    version_list = []

    for gen_coord in gen_coord_list:
        chr_list.append(gen_coord.chromosome)
        start_list.append(gen_coord.start)
        end_list.append(gen_coord.end)
        version_list.append(gen_coord.version)

    return pd.DataFrame(
        {
            'chromosome': chr_list,
            'start': start_list,
            'end': end_list,
            'version': version_list
        }
    )

def from_PairedGenomicCoord_to_DataFrame(paired_gen_coord):
    matched_df = from_GenomicCoord_list_to_DataFrame(
        paired_gen_coord.reference)
    
    matched_df.rename(
        columns={
            'chromosome': 'ref_chromosome',
            'start': 'ref_start',
            'end': 'ref_end',
            'version': 'ref_version'
        }, 
        inplace=True)
    
    matched_df['query_chromosome'] = paired_gen_coord.query.chromosome
    matched_df['query_start'] = paired_gen_coord.query.start
    matched_df['query_end'] = paired_gen_coord.query.end
    matched_df['query_version'] = paired_gen_coord.query.version

    columns_order = [
        'query_chromosome', 'query_start', 'query_end', 'query_version',
        'ref_chromosome', 'ref_start', 'ref_end', 'ref_version'
    ]

    return matched_df.loc[:, columns_order]

def concat_list_of_PairedGenomicCoord_to_DataFrame(paired_gen_coord_list):
    concat_list = []

    for paired_gen_coord in paired_gen_coord_list:
        tmp_df = from_PairedGenomicCoord_to_DataFrame(paired_gen_coord)
        concat_list.append(tmp_df)

    return pd.concat(concat_list).reset_index()
