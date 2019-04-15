#!/usr/bin/env python3

"""
This script mainly holds raw files related methods
"""

import sys
import urllib
import urllib.request
from util import Util


class RawFiles:
    """ This class handles PRIDE submission raw files"""

    base_url = "https://www.ebi.ac.uk/pride/ws/archive/v2/"

    def __init__(self):
        pass

    def download_raw_files_from_ftp(self, accession):
        """
        This method will download all the raw files from PRIDE FTP
        :param accession: PRIDE accession
        :return: None
        """
        request_url = self.base_url + "files/byProject?accession=" + accession + ",fileCategory.value==RAW"
        headers={"Accept": "application/JSON"}

        response = Util.call_api(request_url, headers)

        response_body = response.json()

        for raw_file in response_body:
            ftp_filepath = raw_file['publicFileLocations'][0]['value']
            public_filepath_part = ftp_filepath.rsplit('/', 1)
            print(raw_file['accession'] + " -> " + ftp_filepath)
            new_file_path = raw_file['accession'] + "-" + public_filepath_part[1]
            urllib.request.urlretrieve(ftp_filepath, new_file_path)

def main():
    raw_files =  RawFiles()
    accession = sys.argv[1]
    raw_files.download_raw_files_from_ftp(accession)

if __name__ == '__main__':
    main()