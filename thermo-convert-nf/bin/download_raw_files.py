#!/usr/bin/env python3

"""
This script mainly holds raw files related methods
"""

import requests
import sys
import os
import urllib
import urllib.request

class RawFiles:
    """ This class handles PRIDE submission raw files"""

    base_url = "https://www.ebi.ac.uk/pride/ws/archive/v2/"

    def __init__(self):
        pass

    def downloadRawFiles(self, accession):
        """
        This method will download all the raw files from PRIDE FTP
        :param accession: PRIDE accession
        :return: None
        """
        requestURL = self.base_url + "files/byProject?accession=" + accession + ",fileCategory.value==RAW"
        response = requests.get(requestURL, headers={"Accept": "application/JSON"})

        if (not response.ok) or response.status_code != 200:
            response.raise_for_status()
            sys.exit()

        responseBody = response.json()

        for raw_file in responseBody[:2]:
            ftp_filepath = raw_file['publicFileLocations'][0]['value']
            public_filepath_part = ftp_filepath.rsplit('/', 1)
            print(raw_file['accession'] + " -> " + ftp_filepath)
            new_file_path = raw_file['accession'] + "-" + public_filepath_part[1]
            urllib.request.urlretrieve(ftp_filepath, new_file_path)

def main():
    raw_files =  RawFiles()
    accession = sys.argv[1]
    raw_files.downloadRawFiles(accession)

if __name__ == '__main__':
    main()