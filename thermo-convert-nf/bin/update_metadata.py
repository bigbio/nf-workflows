#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script mainly fetches a access token from the API and
secondly update the msrun metadata via API endpoint
"""

import sys
import json
import os
import urllib
import requests
from authentication import Authentication
from filehandling import FileHandling

class MsRun:
    base_url = "https://www.ebi.ac.uk/pride/ws/archive/v2/"

    def __init__(self):
        pass

    def updateMsrunMetadata(self, filename, token):
        # get project file accession from the prefix of the file name (e.g: PXF00000145820)
        accession = filename.split('-', 1)[0]

        url = self.base_url + "msruns/" + accession + "/updateMetadata"
        headers = {'Content-type': 'application/json', 'Accept': 'application/json', 'Authorization': 'Bearer ' + token}

        with open(filename) as json_file:
            data = json.load(json_file)
            print(json.dumps(data))
            response = requests.put(url, data=json.dumps(data), headers=headers)

            if (not response.ok) or response.status_code != 200:
                response.raise_for_status()
                sys.exit()
            else:
                print(response)

def main():
    filename = sys.argv[1]
    username = sys.argv[2]
    password = sys.argv[3]

    # Get user token to make calls with PRIDE API
    authentication=Authentication()
    token = authentication.get_token(username, password)

    # Format extracted metadata to compatible with PRIDE API endpoint
    fileHandling=FileHandling()
    fileHandling.wrap_with_ms_run_metadata(filename)

    # Update msrun metatdata
    msrun = MsRun()
    msrun.updateMsrunMetadata(filename, token)

if __name__ == '__main__':
    main()
