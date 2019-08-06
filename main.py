from convert_annotation import read
from convert_annotation import formatter

def main(map_csv_path, query_version, ref_version, 
    query_path, out_csv_path, description=''):
    cc = read.from_CSV_to_ConvertCoordinates(
        map_csv_path, query_version, ref_version, description='')

    querys = read.parse_query_list(query_path)

    paired_gen_coord_list = cc.recursively_get_dmel_coordinates(
        query_version, ref_version, querys)

    out_df = formatter.concat_list_of_PairedGenomicRanges_to_DataFrame(
        paired_gen_coord_list)

    out_df.to_csv(out_csv_path)
