#!/usr/bin/env python
# coding: utf-8

import requests

## set up to automatically upload to existing zenodo
#deposit = '6193041'

ACCESS_TOKEN = 'SNMWz9fad8X87I4JVL51D664P5QFWUCNYHDMPbNsiEnql0OD76pXsDm4yn4C'

### get existing
#params = {'access_token': ACCESS_TOKEN}
#r = requests.get('https://zenodo.org/api/deposit/depositions/' + deposit,
#                   params=params,
#                   json={},
                   # Headers are not necessary here since "requests" automatically
                   # adds "Content-Type: application/json", because we're using
                   # the "json=" keyword argument
                   # headers=headers,
                   # headers=headers)
#                 )

### create new
params = {'access_token': ACCESS_TOKEN}
r = requests.post('https://zenodo.org/api/deposit/depositions',
                   params=params,
                   json={},
                   # Headers are not necessary here since "requests" automatically
                   # adds "Content-Type: application/json", because we're using
                   # the "json=" keyword argument
                   # headers=headers,
                   # headers=headers)
                 )

print(r.status_code)

bucket_url = r.json()["links"]["bucket"]

filename = "restart.tar.gz"
path = "/glade/campaign/univ/uwas0070/vcooper/waveice/zenodo_modeloutput/minfiles/%s" % filename

with open(path, "rb") as fp:
    r = requests.put(
        "%s/%s" % (bucket_url, filename),
        data=fp,
        params=params,
    )

print(r.status_code)
    
