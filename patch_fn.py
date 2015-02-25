#!/usr/bin/env python

import dxencode
import dxpy
import json
import requests

def patch(obs):
    for fob in obs:
        if fob['file_format'] == 'fastq' or fob['status'] == 'revoked':
            continue
        fn = fob['submitted_file_name']
        folder = dxpy.describe(dxpy.find_one_data_object(name=fn.strip('/'), project='project-BQkYKg00F1GP55qQ9Qy00VP0')['id'])['folder']
        newfn = folder+'/'+fn.strip('/')
        print "Patch: %s with %s" % (fn, newfn)
        res = requests.patch(srv+fob['@id'], auth=(id,pw), data=json.dumps({'submitted_file_name': newfn}),headers={'content-type': 'application/json'})
        try:
            res.raise_for_status()
            print "Success"
        except Exception, e:
            print "Failed %s" % e


(id,pw,srv) = dxencode.processkey('www')

accs = dxencode.encoded_get(srv+'ENCSR000AEV', AUTHID=id, AUTHPW=pw).json()['original_files']
file_objs = [ dxencode.encoded_get(srv+acc, AUTHID=id, AUTHPW=pw).json() for acc in accs ]

patch(file_objs)