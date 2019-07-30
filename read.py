import pandas as pd
from convert_annotation.classes import ConvertCoordinates

def from_CSV_to_ConvertCoordinates(csv_path, version1, version2, description=''):
    """ Read CSV file into ConvertCoordinates object. """
    return ConvertCoordinates(
        pd.read_csv(csv_path), version1, version2, description)

