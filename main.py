from convert_annotation import read
from convert_annotation import formatter

def main(map_csv_path, version1, version2, query_path, out_csv_path, description=''):
    cc = read.from_CSV_to_ConvertCoordinates(
        map_csv_path, version1, version2, description='')

    querys = read.parse_query_list(query_path)

    gen_coord_list = cc.recursive_find_dmel_coordinates(querys)

    out_df = formatter.from_GenomicCoord_list_to_DataFrame(gen_coord_list)


