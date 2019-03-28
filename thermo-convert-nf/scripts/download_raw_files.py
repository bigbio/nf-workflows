#!/usr/bin/env python
import requests
import sys
import os
import urllib

accession = sys.argv[1]
requestURL = "https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession=" + accession + ",fileCategory.value==RAW"

response = requests.get(requestURL, headers={"Accept": "application/JSON"})

if (not response.ok) or response.status_code != 200:
    response.raise_for_status()
    sys.exit()

responseBody = response.json()

for raw_file in responseBody[:2]:
    ftp_filepath = raw_file['publicFileLocations'][0]['value']
    public_filepath_part = ftp_filepath.rsplit('/', 1)
    print(raw_file['accession'] + " -> "+ ftp_filepath)
    new_file_path=raw_file['accession'] + "-"+  public_filepath_part[1]
    urllib.urlretrieve(ftp_filepath, new_file_path)