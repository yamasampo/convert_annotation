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
