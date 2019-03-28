#!/usr/bin/env python
import requests
import sys
import os

accession = sys.argv[1]
requestURL = "https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession=" + accession + ",fileCategory.value==RAW"

response = requests.get(requestURL, headers={"Accept": "application/JSON"})

if (not response.ok) or response.status_code != 200:
    response.raise_for_status()
    sys.exit()

responseBody = response.json()
raw_file = responseBody[0]
public_filepath = raw_file['publicFileLocations'][0]['value']
public_filepath_part = public_filepath.rsplit('/', 1)
public_filepath_folder = public_filepath_part[0] + "/"
text_file = open("project_path.txt", "w")
text_file.write(public_filepath_folder)
text_file.close()