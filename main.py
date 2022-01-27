from . import read
from . import formatter

def main(
        map_csv_path:str, 
        query_version:str, 
        ref_version:str, 
        query_path:str, 
        out_csv_path:str, 
        **kwargs
        ):
    """ Outputs a file with converted coordinates from a given query file. 
    
    Parameters
    ----------
    map_csv_path: str
        A path to a file including the corresponding coordinates. 
    query_version: str, int or tuple
        A key to keep track of query data.
    ref_version: str, int or tuple
        A key to keep track of reference data. 
    query_path: str
        A path to a file including that a user wishes to convert its cooredinate
        to the reference version. 
    out_csv_path: str
        A path to an output file. 
    
    """
    cc = read.from_CSV_to_ConvertCoordinates(
        map_csv_path, query_version, ref_version)

    querys = read.parse_query_list(query_path, **kwargs)

    paired_gen_coord_list = cc.recursively_get_dmel_coordinates(
        query_version, ref_version, querys)

    out_df = formatter.concat_list_of_PairedGenomicRanges_to_DataFrame(
        paired_gen_coord_list)

    out_df.to_csv(out_csv_path)
