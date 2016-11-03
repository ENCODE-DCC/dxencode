import os, sys, json
import subprocess, commands, requests, urlparse
from datetime import datetime
import time
#import shlex

import logging

ENCODED_VERSION = "1"

INTERNAL_STATUS_BLOCKS = ["requires lab review", "unrunnable"]
'''Experients with these internal_statuses should not be assembled or launched.'''

KEYFILE = 'keypairs.json'  ## see processkey() Note this file must be in gitignore!
PRODUCTION_SERVER = 'https://www.encodeproject.org'
DEFAULT_SERVER = PRODUCTION_SERVER
S3_SERVER='s3://encode-files/'

DCC_PIPELINE_LAB_NAME='encode-processing-pipeline'
DCC_PIPELINE_LAB='/labs/'+DCC_PIPELINE_LAB_NAME+'/'
DEFAULT_DCC_AWARD='/awards/U41HG006992/'

logger = logging.getLogger('encd') # Callers should either use dxencode.logger or set dxencode.logger = local.logger()

SAVED_KEYS = {}

# TODO: make Encd object and replace prime_server_key with self.server_key
prime_server_key = "default" # ordinarily is set to the first non-default key seen.

def set_server_key(key):
    '''Explicitly sets prime_server_key, which ordinarily is set to the first non-default key seen.'''
    global prime_server_key
    prime_server_key = key

def get_server_key():
    '''Returns the server_key.'''
    global prime_server_key
    return prime_server_key

def get_server():
    '''Returns the server.'''
    (authid,authpw,server) = find_keys()
    return server

def find_keys(server=None,authid=None,authpw=None):
    '''returns authentication keys.'''
    if server == None or authid == None or authpw == None:
        return processkey( server ) # Assume server is server_key or None to get default
    return (authid,authpw,server)   # pass through

    

def processkey(key):
    ''' check encodedD access keys; assuming the format:
    {
        "default":
                {"server":"https://www.encodeproject.org", "key":"rand_user_name", "secret":"rand_password"},
        "test":
                {"server":"https://test.encodedcc.org", "key":"rand_user_name", "secret":"rand_password"},
        "www":
                {"server":"https://www.encodeproject.org", "key":"rand_user_name", "secret":"rand_password"}
    }
    '''
    global prime_server_key
    if key == None:
        key = prime_server_key
    if key in SAVED_KEYS:
        return SAVED_KEYS[key]

    if key:
        keysf = open(KEYFILE,'r')
        keys_json_string = keysf.read()
        keysf.close()
        keys = json.loads(keys_json_string)
        key_dict = keys[key]
    else:
        key_dict = {}
    AUTHID = key_dict.get('key')
    AUTHPW = key_dict.get('secret')
    if key:
        SERVER = key_dict.get('server')
    else:
        SERVER = 'https://www.encodeproject.org/'

    if not SERVER.endswith("/"):
        SERVER += "/"

    SAVED_KEYS[key] = (AUTHID,AUTHPW,SERVER)
    # Gets set to first non-default key
    if prime_server_key == "default":
        prime_server_key = key
    return (AUTHID,AUTHPW,SERVER)
    ## TODO possibly this should return a dict
    
def post_obj(obj_type,obj_meta, SERVER=None, AUTHID=None, AUTHPW=None):
    ''' Posts a json object of a given type to the encoded database. '''
    (AUTHID,AUTHPW,SERVER) = find_keys(SERVER, AUTHID, AUTHPW)
    HEADERS = {
        'Content-type': 'application/json',
        'Accept': 'application/json',
    }
    r = requests.post(
        SERVER + obj_type,
        auth=(AUTHID, AUTHPW),
        data=json.dumps(obj_meta),
        headers=HEADERS,
    )
    try:
        r.raise_for_status()
    except:
        logger.error('Submission failed: %s %s' % (r.status_code, r.reason))
        logger.error(r.text)
        raise
    item = r.json()['@graph'][0]
    # TODO: Look for returned["status"] == "success"
    #returned = r.json()
    #if "status" not in returned or returned["status"] != "success":
    #    print >> sys.stderr, "* WARNING: request to post '%s'..." % obj_type
    #    print >> sys.stderr, "Posting " + SERVER + obj_type
    #    print >> sys.stderr, json.dumps(obj_meta, indent=4, sort_keys=True)
    #    print >> sys.stderr, "Returned:"
    #    print >> sys.stderr, json.dumps(item, indent=4, sort_keys=True)
    ##    print >> sys.stderr, json.dumps(r.json(), indent=4, sort_keys=True)
    return item

def patch_obj(obj_id, obj_meta, SERVER=None, AUTHID=None, AUTHPW=None):
    ''' Patches a json object of a given type to the encoded database. '''
    (AUTHID,AUTHPW,SERVER) = find_keys(SERVER, AUTHID, AUTHPW)
    #HEADERS = { 'Content-type': 'application/json' }
    HEADERS = {
        'Content-type': 'application/json',
        'Accept': 'application/json',
    }
    r = requests.patch(
        SERVER + obj_id,
        auth=(AUTHID, AUTHPW),
        data=json.dumps(obj_meta),
        headers=HEADERS,
    )
    try:
        r.raise_for_status()
    except:
	    #if not r.status_code == 200:
		#    print >> sys.stderr, r.text
        logger.error('Patch of %s failed: %s %s' % (obj_id, r.status_code, r.reason))
        logger.error(r.text)
        raise

    item = r.json()['@graph'][0]
    #print >> sys.stderr, "* request to patch %s to %s..." % (obj_id,SERVER)
    #print >> sys.stderr, json.dumps(item, indent=4, sort_keys=True)
    return item
           

def post_file(filename, file_meta, SERVER=None, AUTHID=None, AUTHPW=None):
    ''' take a file object on local file system, post meta data and cp to AWS '''
    (AUTHID,AUTHPW,SERVER) = find_keys(SERVER, AUTHID, AUTHPW)
    HEADERS = {
        'Content-type': 'application/json',
        'Accept': 'application/json',
    }

    # Handle case where file_obj is already there but file is not!
    if file_meta.get('accession') != None:
        logger.debug("Patching file object.")
        assert file_meta['status'] == "upload failed"
        file_meta['status'] = "uploading"  # Laurence will change encodeD function to take care of this
        item = patch_obj('files/'+file_meta.get('accession'),file_meta, SERVER, AUTHID, AUTHPW)
        item = post_obj('files/'+file_meta.get('accession')+"/@@upload",{}, SERVER, AUTHID, AUTHPW)
    else:
        logger.debug("Posting file object.")
        assert file_meta.get('status') == None
        item = post_obj('file',file_meta, SERVER, AUTHID, AUTHPW)

    # Could look for returned["status"] == "success"
    #returned = r.json()
    #if "status" not in returned or returned["status"] != "success" or \
    if 'upload_credentials' not in item:
        print >> sys.stderr, "* ERROR: request to post %s to %s..." % (filename,SERVER)
        print >> sys.stderr, json.dumps(file_meta, indent=4, sort_keys=True)
        print >> sys.stderr, "* Returned..."
        print >> sys.stderr, json.dumps(item, indent=4, sort_keys=True)

    # if cred missing, just let it fail on the next statement
    creds = item['upload_credentials']
    env = os.environ.copy()
    env.update({
        'AWS_ACCESS_KEY_ID': creds['access_key'],
        'AWS_SECRET_ACCESS_KEY': creds['secret_key'],
        'AWS_SECURITY_TOKEN': creds['session_token'],
    })

    # POST file to S3
    logger.debug("Uploading file.")
    start = datetime.now()
    try:
        subprocess.check_call(['aws', 's3', 'cp', filename, creds['upload_url']], env=env)
    except:
        logger.debug("Retry 2/3 uploading in 1 minute...")
        time.sleep(60)
        try:
            subprocess.check_call(['aws', 's3', 'cp', filename, creds['upload_url']], env=env)
        except:
            logger.debug("Retry 3/3 uploading in 3 minutes...")
            time.sleep(180)
            try:
                subprocess.check_call(['aws', 's3', 'cp', filename, creds['upload_url']], env=env)
            except:
                logger.debug("Upload failed")
                # Try to set status to "upload failed"
                file_meta['status'] = "upload failed"
                item = patch_obj('files/'+item.get('accession'),file_meta, SERVER, AUTHID, AUTHPW)
                assert item['status'] == "upload failed"
                #raise Don't raise exception on half error... the accession needs to be added to the dx file.

    if item.get('status','uploading') != "upload failed":
        #assert item['status'] == "uploading"
        end = datetime.now()
        duration = end - start
        logger.debug("Uploaded in %.2f seconds" % duration.seconds)
            
    return item


def get_object(url, AUTHID=None, AUTHPW=None):
    ''' executes GET on Encoded server without without authz '''
    ##TODO possibly add try/except looking for non 4xx?
    HEADERS = {'content-type': 'application/json'}
    if AUTHID and AUTHPW:
        response = requests.get(url, auth=(AUTHID,AUTHPW), headers=HEADERS)
    else:
        response = requests.get(url, headers=HEADERS)
    return response


def lookup_json(path, key=None, frame='object', must_find=False):
    '''Commonly used method to get a json object from encodeD.'''
    (AUTHID,AUTHPW,SERVER) = find_keys(key)
    url = SERVER + path + '/?format=json&frame=' + frame
    #print >> sys.stderr, url
    response = get_object(url, AUTHID, AUTHPW)
    try:
        response.raise_for_status()
        json_obj = response.json()
    except:
        if must_find:
            print >> sys.stderr, "Path to json object '%s' not found." % path
            print >> sys.stderr, 'Lookup failed: %s %s' % (response.status_code, response.reason)
            sys.exit(1)
        return None
    return json_obj


def get_bucket(f_obj, SERVER=None, AUTHID=None, AUTHPW=None):
    ''' returns aws s3 bucket and file name from encodeD file object (f_obj)'''
    #make the URL that will get redirected - get it from the file object's href property
    (AUTHID,AUTHPW,SERVER) = find_keys(SERVER, AUTHID, AUTHPW)
    encode_url = urlparse.urljoin(SERVER,f_obj.get('href'))
    logger.debug(encode_url)

    #stream=True avoids actually downloading the file, but it evaluates the redirection
    r = requests.get(encode_url, auth=(AUTHID,AUTHPW), headers={'content-type': 'application/json'}, allow_redirects=True, stream=True)
    try:
        r.raise_for_status
    except:
        logger.error('%s href does not resolve' %(f_obj.get('accession')))
        sys.exit()
    logger.debug(r)

    #this is the actual S3 https URL after redirection
    s3_url = r.url
    logger.debug(s3_url)

    #release the connection
    r.close()

    #split up the url into components
    o = urlparse.urlparse(s3_url)
    opath = o.path.replace("/http://encode-files.s3.amazonaws.com", "") # Hacked to make sure Aditi's assemble works

    #pull out the filename
    filename = os.path.basename(o.path)

    #hack together the s3 cp url (with the s3 method instead of https)
    return filename, S3_SERVER.rstrip('/') + opath

def file_in_list(looking_for_file,file_list):
    md5 = looking_for_file.get('md5sum')
    for a_file in file_list:
        if md5 == a_file.get('md5sum'):
            return True
    return False

def files_to_map(exp_obj):
    if not exp_obj or not exp_obj.get('files'):
        return []
    else:
        files = []
        for file_obj in exp_obj.get('files'):
            if (file_obj.get('output_type') == 'reads' or file_obj.get('output_type') == 'raw data') and \
               file_obj.get('file_format') == 'fastq' and \
               file_obj.get('replicate') and file_obj.get('replicate').get('biological_replicate_number') and \
                                             file_obj.get('replicate').get('technical_replicate_number') and \
               not file_in_list(file_obj,files):
               files.extend([file_obj])
            elif file_in_list(file_obj,files):
                logger.warning('%s:%s Duplicate file md5sum, ignoring.' %(exp_obj.get('accession'),file_obj.get('accession')))
                print >> sys.stderr, 'WARNING: %s:%s Duplicate filename, ignoring.' % \
                                                                        (exp_obj.get('accession'),file_obj.get('accession'))
                #return []
        return files

def replicates_to_map(experiment, files):
    if not files:
        return []
    else:
        reps_with_files = set([ f['replicate']['uuid'] for f in files if f.get('replicate') ])
        return [ r for r in experiment['replicates'] if r['uuid'] in reps_with_files ]

def is_paired_ended(experiment):
    ''' this is likely not the most efficient way to do this'''

    mapping = choose_mapping_for_experiment(experiment, warn=True)
    reps_paired = {}
    for rep in mapping.keys():
        p = mapping[rep].get('paired', [])
        up = mapping[rep].get('unpaired', [])
        # I should be a CS guy and do some XOR thing here.
        if p and not up:
            reps_paired[rep] = True
        elif up and not p:
            reps_paired[rep] = False
        else:
            print >> sys.stderr, "Mixed mapping for replicate %s/%s" %(experiment['accession'], rep)
            print >> sys.stderr, "Paired: %s" % ([ f['accession'] for f in p ])
            print >> sys.stderr, "Unpaired: %s" % ([ f['accession'] for f in up ])
            sys.exit(1)

    trues = len([ v for v in reps_paired.values() if v ])
    falses = len([ v for v in reps_paired.values() if not v])
    if trues and falses:
        print >> sys.stderr, "Mixed mapping for replicates in experiment %s" % (experiment['accession'])
        print >> sys.stderr, reps_paired
    else:
        if trues:
            return True
        elif falses:
            return False

    print >> sys.stderr, "Never get here"
    sys.exit(1)

def choose_mapping_for_experiment(experiment,warn=True):
    ''' for a given experiment object, fully embedded, return experimental info needed for mapping
        returns an dict keyed by [biological_rep][technical_rep]
        with information for mapping (sex, organism, paired/unpaired files, library id)
    '''

    exp_id = experiment['accession']
    exp_files = files_to_map(experiment)
    files = [f for f in exp_files if f.get('replicate') and 
                                     f.get('replicate').get('biological_replicate_number') and
                                     f.get('replicate').get('technical_replicate_number')]
    replicates = replicates_to_map(experiment, files)
    mapping = {}

    if files:
        for rep in replicates:
            biorep_n = rep.get('biological_replicate_number')
            techrep_n = rep.get('technical_replicate_number')


            try:
                library = rep['library']['accession']
                sex = rep['library']['biosample'].get('sex', "male")
                if sex != "male" and sex != "female":
                    if warn:
                        print >> sys.stderr, "WARN: using male replacement for %s" % sex
                    sex = "male"
                organism = rep['library']['biosample']['donor']['organism']['name']
            except KeyError:
                print >> sys.stderr, "Error, experiment %s replicate %s_%s missing info\n%s" % (exp_id,biorep_n,techrep_n,rep)
                sys.exit(0)

            rep_files = [f for f in files if f.get('replicate').get('biological_replicate_number') == biorep_n and
                                                f.get('replicate').get('technical_replicate_number') == techrep_n ]
            paired_files = []
            unpaired_files = []
            while rep_files:
                file_object = rep_files.pop()
                if file_object.get('paired_end') == None: # group all the unpaired reads for this biorep together
                    unpaired_files.extend([ file_object ])
                elif file_object.get('paired_end') in ['1','2']:
                    if file_object.get('paired_with'):
                        mate = next((f for f in rep_files if f.get('@id') == file_object.get('paired_with')), None)
                    else: #have to find the file that is paired with this one
                        mate = next((f for f in rep_files if f.get('paired_with') == file_object.get('@id')), None)
                    if mate:
                        rep_files.remove(mate)
                    elif warn:
                        logger.warning('%s:%s could not find mate' %(experiment.get('accession'), file_object.get('accession')))
                        mate = {}
                    paired_files.extend([ (file_object, mate) ])

            mapping[(biorep_n, techrep_n)] = {
                "library": library,
                "sex": sex,
                "organism": organism,
                "paired": paired_files,
                "unpaired": unpaired_files,
                "replicate_id": rep['@id']
            }
            if rep_files and warn:
                logger.warning('%s: leftover file(s) %s' % (exp_id, rep_files))
    elif warn:
        logger.warning('%s: No files to map' % exp_id)
    return mapping

def get_exp(experiment,must_find=True,warn=False,key=None):
    '''Returns all replicate mappings for an experiment from encoded.'''

    (AUTHID,AUTHPW,SERVER) = find_keys(key)
    url = SERVER + 'experiments/%s/?format=json&frame=embedded' % experiment
    try:
        response = get_object(url, AUTHID, AUTHPW)
        exp = response.json()
    except:
        if must_find:
            print >> sys.stderr, "Experiment %s not found." % experiment
            print >> sys.stderr, "Auth %s %s may not be valid for %s." % (AUTHID, AUTHPW, url)
            #print >> sys.stderr, response #.json()
            sys.exit(1)
        return None
    if exp == None or exp["status"] == "error":
        if must_find:
            print >> sys.stderr, "Experiment %s not found." % experiment
            print >> sys.stderr, response.json()
            sys.exit(1)
        return None

    return exp

def get_assay_type(experiment,exp=None,key=None,must_find=True,warn=False):
    '''Looks up encoded experiment's assay_type, normalized to lower case.'''
    if exp == None:
        exp = get_exp(experiment,key=key,must_find=must_find,warn=warn)
    if "assay_term_name" not in exp:
        if must_find:
            print >> sys.stderr, "No 'assay_term_name' found for experiment %s." % experiment
            sys.exit(1)
        return None
    if exp["assay_term_name"] in [  "RNA-seq", \
                                    "shRNA knockdown followed by RNA-seq", \
                                    "CRISPR genome editing followed by RNA-seq", \
                                    "single cell isolation followed by RNA-seq", \
                                    "siRNA knockdown followed by RNA-seq" ]:
        #if exp["replicates"][0]["library"]["size_range"] in [">200", "300-350", "350-450"]:
        # Now more: "150-400","149-512","151-499","153-499","157-497"        
        size_range = exp["replicates"][0]["library"]["size_range"]
        if size_range.startswith('>'):
            try:
                min_size = int(size_range[1:])
            except:
                min_size = 0
            max_size = min_size
        else:
            try:
                sizes = size_range.split('-')
                min_size = int(sizes[0])
                max_size = int(sizes[1])
            except:
                min_size = 0
                max_size = 0
        if max_size <= 200 and max_size != min_size:
            return "small-rna-seq"
        elif min_size >= 150:
            return "long-rna-seq"
        elif (min_size + max_size)/2 >= 235: # This is some wicked voodoo (SRNA:108-347=227; LRNA:155-315=235)        
            return "long-rna-seq"
        elif min_size == 120 and max_size == 200: # Another ugly exception!        
            return "long-rna-seq"
        else:
            return "small-rna-seq"
    elif exp["assay_term_name"] == "whole-genome shotgun bisulfite sequencing" \
      or exp["assay_term_name"] == "shotgun bisulfite-seq assay":
        return "dna-me"
    elif exp["assay_term_name"] == "CAGE":
        return "rampage"
    #elif exp["assay_term_name"] == "RAMPAGE":
    #    return "rampage"
    #elif exp["assay_term_name"] == "ChIP-seq":
    #    return "chip-seq"

    return exp["assay_term_name"].lower()

def get_exp_type(experiment,exp=None,supported_types=None):
    '''Looks up encoded experiment's assay_type, normalized to known supported tokens.'''
    if exp == None:
        exp = get_exp(experiment)
    exp_type = get_assay_type(experiment,exp)

    if supported_types != None and exp_type not in supported_types:
        print >> sys.stderr, "Experiment %s has unsupported assay type of '%s'" % (experiment,exp_type)
        return None
    return exp_type

def get_full_mapping(experiment,exp=None,key=None,must_find=True,warn=False):
    '''Returns all replicate mappings for an experiment from encoded.'''

    if exp == None:
        exp = get_exp(experiment,key=key,must_find=must_find,warn=warn)

    if not exp.get('replicates') or len(exp['replicates']) < 1:
        if must_find:
            print >> sys.stderr, "No replicates found in %s" % experiment
            sys.exit(1)
        return None

    return choose_mapping_for_experiment(exp,warn=warn)

def get_replicate_mapping(experiment,biorep=None,techrep=None,full_mapping=None,key=None, \
                                                                                    must_find=True):
    '''Returns replicate mappings for an experiment or specific replicate from encoded.'''
    if full_mapping == None:
        full_mapping = get_full_mapping(experiment,key=key,must_find=must_find,warn=False)

    try:
        return full_mapping[(biorep,techrep)]
    except KeyError:
        if must_find:
            print >> sys.stderr, "Specified replicate: rep%s_%s could not be found in mapping of %s." % \
                                                                                            ( biorep, techrep, experiment )
            print >> sys.stderr, json.dumps(full_mapping,indent=4)
            sys.exit(1)
    return None

def get_reps(exp_id, load_reads=False, exp=None, full_mapping=None, key=None):
    '''For a given exp_id (accession) returns a "rep" list as used by assemble, launch, etc.'''
        
    reps = []
    # Must look through exp and find all replicates!
    if full_mapping == None:
        full_mapping = get_full_mapping(exp_id,exp,key=key)
    if full_mapping != None:
        for (br,tr) in sorted( full_mapping.keys() ):
            rep = { 'br': br, 'tr': tr,'rep_tech': 'rep' + str(br) + '_' + str(tr) }
            mapping = get_replicate_mapping(exp_id,br,tr,full_mapping) # must_find
            rep['organism'] = mapping['organism']
            rep['sex'] = mapping['sex']
            rep['library_id'] = mapping['library']
            rep['replicate_id'] = mapping['replicate_id']
            rep['paired_end'] = False  # Default in case read files are missing
            rep['has_reads'] = False
            if mapping['paired'] and not mapping['unpaired']:
                rep['paired_end'] = True
                rep['has_reads'] = True
            elif mapping['unpaired'] and not mapping['paired']:
                rep['paired_end'] = False
                rep['has_reads'] = True
            elif mapping['paired'] and mapping['unpaired']:
                print >> sys.stderr, "ERROR: Replicate %d_%d has both paired(%s) and unpaired(%s) reads." % \
                    (br,tr,len(mapping['paired']), len(mapping['unpaired']))
                #print >> sys.stderr, json.dumps(mapping,indent=4)
                sys.exit(1)                
                rep['paired_end'] = False
                rep['has_reads'] = True
                
            # Load read and control files, only if requested.
            run_type = None  # The files should be consistent as all single-end or paired-end, but some mapping got it wrong
            if load_reads and rep['has_reads']:
                rep['fastqs'] = { "1": [], "2": [] }
                rep['controls'] = []
                if rep['paired_end']:
                    for (p1, p2) in mapping['paired']:
                        if p1['status'] not in ['released','in progress']:
                            continue
                        rep['fastqs'][p1['paired_end']].append(p1['accession']+".fastq.gz")
                        if "run_type" in p1:
                            if run_type == None or run_type == "single-ended":
                                run_type = p1["run_type"]
                            #elif run_type != p1["run_type"]:
                            #    # Warn?
                        if p2 != None and 'paired_end' in p2:
                            rep['fastqs'][p2['paired_end']].append(p2['accession']+".fastq.gz")
                        if 'controlled_by' in p1:
                            rep['controls'].append( p1['controlled_by'] )
                        if p2 != None and 'controlled_by' in p2:
                            rep['controls'].append( p2['controlled_by'] )
                else: # not rep['paired_end']:
                    for f in mapping['unpaired']:
                        if f['status'] not in ['released','in progress']:
                            continue
                        rep['fastqs']['1'].append( f['accession']+".fastq.gz" )
                        if "run_type" in f:
                            if run_type == None or run_type == "single-ended":
                                run_type = f["run_type"]
                    if 'controlled_by' in mapping['unpaired']:
                        rep['controls'].append( mapping['unpaired']['controlled_by'] )
                    elif 'controlled_by' in mapping['unpaired'][0]:
                        rep['controls'].append( mapping['unpaired'][0]['controlled_by'] )
                if len(rep['controls']) == 0:
                    rep['controls'] = None
            
            # One more test because of non-standard data: single-end data reporting 'paired' == 1 !
            if rep['paired_end'] and run_type == "single-ended" and len(rep['fastqs']['2']) == 0:
                rep['paired_end'] = False
                
            reps.append( rep )

    return reps

def get_file(file_acc,must_find=False,key=None):
    '''Returns file object from encoded, given an accession.'''

    (AUTHID,AUTHPW,SERVER) = find_keys(key)

    if file_acc.startswith("/files/") and file_acc.endswith('/'):
        url = SERVER + '%s?format=json&frame=embedded' % file_acc
    else:
        url = SERVER + '/files/%s/?format=json&frame=embedded' % file_acc
    try:
        response = get_object(url, AUTHID, AUTHPW)
        file_obj = response.json()
    except:
        if must_find:
            print >> sys.stderr, "File %s not found." % file_acc
            sys.exit(1)
        return None

    return file_obj

def get_exp_files(exp_obj,output_types=[],lab=None,key=None):
    '''Returns list of file objs associated with an experiment, filtered by zero or more output_types.'''
    if not exp_obj or not exp_obj.get('files'):
        return []
    files = []
    accessions = []
    #print >> sys.stderr, "DEBUG: Found %d original_files" % len(exp_obj['original_files']) 
    for file_acc in exp_obj['original_files']:
        file_obj = None
        acc = file_acc[7:18]
        if acc in accessions:
            continue
        accessions.append(acc)
        for f_obj in exp_obj['files']:
            if acc == f_obj['accession']:
                file_obj = f_obj
                break
        if file_obj == None:
            file_obj = get_file(file_acc,key=key)
        if file_obj == None:
            continue
        #print >> sys.stderr, " * Found: %s [%s] status:%s %s" % \
        #                  (file_obj['accession'],file_obj['output_type'],file_obj['status'],file_obj['submitted_file_name'])
        out_type = file_obj.get('output_type')
        #print >> sys.stderr, "DEBUG: File %s is of out_type '%s'" % (acc,out_type) 
        if len(output_types) > 0 and (out_type == None or out_type not in output_types):
            continue
        if lab != None:
            file_lab = file_obj.get('lab')
            if file_lab != None:
                if (isinstance(file_lab,unicode) or isinstance(file_lab,str)) and file_lab != "/labs/"+lab+"/":
                    continue
                elif "name" in file_lab and file_lab["name"] != lab:
                    #print >> sys.stderr, "DEBUG: File %s is of out_type '%s' and lab '%s' not '%s'" % (acc,out_type,file_lab["name"],lab) 
                    continue
        if file_obj.get('status') not in ["released","uploaded","uploading","in progress"]: # further restricted by caller.
            continue
        if not file_in_list(file_obj,files):
           files.extend([file_obj])
    return files

def exp_is_pe(exp,exp_files=None,rep_tech=None,server_key=None):
    '''Determine if this experiment is expected to be 'paired-end' as opposed to 'single-end'.'''
    
    if rep_tech != None:
        reps = get_reps(exp_id=exp.get('accession'), load_reads=False, exp=exp, full_mapping=None, key=server_key)
        for rep in reps:
            if rep.get('rep_tech') == rep_tech:
                if rep.get('paired_end') == False:
                    return False  # no more is needed!
        
    # But if rep_tech not requested or not found, OR found to be paired, then must check all fastqs
    found_fastq = False
    all_fastqs_pe = True
    
    if exp_files == None:
        exp_files = get_exp_files(exp,key=server_key)
        
    for f_obj in exp_files:
        if f_obj.get("file_format") != "fastq":
            continue
        found_fastq = True
        if f_obj.get("paired_with") == None and f_obj.get("run_type") != "paired-ended":
            all_fastqs_pe = False
            
    return (found_fastq and all_fastqs_pe) 

def select_alias(aliases, prefix='dnanexus:',must_find=True):
    '''returns the alias of a given prefix'''
    for alias in aliases:
        if alias.startswith(prefix):
            return alias
    if must_find:
        print >> sys.stderr, "ERROR: Failed to find alias of prefix: '"+prefix+"' in:"
        print >> sys.stderr, aliases
        sys.exit(1)
    
    return None


def rep_is_umi(exp,rep=None,exp_files=None,rep_tech=None,server_key=None):
    '''Determine if this technical replicate is has fastqs all marked as UMI or all non-UMI.'''
    
    if rep == None:
        reps = get_reps(exp_id=exp.get('accession'), load_reads=False, exp=exp, full_mapping=None, key=server_key)
        for one_rep in reps:
            if one_rep.get('rep_tech') == rep_tech:
                rep = one_rep
                break

    if exp_files == None:
        exp_files = get_exp_files(exp,key=server_key)
    #print >> sys.stderr, ">  Looking through %d files for rep%d_%d %s" % (len(exp_files),rep.get('br',0),rep.get('tr',0),rep.get('replicate_id'))
        
    umi_found_true  = False
    umi_found_false = False
    umi = None
    
    barcodes = []
    for f_obj in exp_files:
        if f_obj.get("file_format") != "fastq":
            continue
        if f_obj["replicate"]['@id'] != rep.get('replicate_id'):
            continue
        #print >> sys.stderr, ">>  Looking for umi on %s.%s %s with %d flow objects" % \
        #    (f_obj.get("accession"),f_obj.get("file_format"),f_obj["replicate"]['@id'],len(f_obj.get("flowcell_details",[])))
        umi_found = False
        for flow in f_obj.get("flowcell_details",[]):
            barcode = flow.get('barcode')
            if barcode == None:
                continue
            barcodes.append(barcode)
            if barcode == "UMI" or barcode.startswith("SSLIB"):
                umi_found = True
                umi_found_true  = True
        if not umi_found :
            umi_found_false = True
        
    if umi_found_true and umi_found_false:
        print >> sys.stderr, "ERROR: mixed UMI setting on the fastqs for replicate %s." % rep.get('replicate_id')
        sys.exit(1)
    if not umi_found_true and not umi_found_false:
        print >> sys.stderr, "WARNING: could not detect UMI for replicate %s." % rep.get('replicate_id')
    #print >> sys.stderr, "Found: UMI for replicate %s to be %s." % (rep.get('replicate_id'),str(umi_found_true))
    return (umi_found_true, barcodes) 


def exp_patch_internal_status(exp_id, internal_status, key=None, test=False):
    '''Updates encodeD Experiment with an internal status.'''
    
    allowed = [ 'pipeline ready', 'processing', 'pipeline completed', 'requires lab review' ]
    
    if internal_status not in allowed:
        print >> sys.stderr, "  * ERROR: Attempting to set internal status of %s to '%s'." % (exp_id,internal_status)
        return False
    if not exp_id.startswith('ENCSR'):
        print >> sys.stderr, "  * ERROR: Attempting to set internal status on experiment without valid accession %s." % exp_id
        return False
        
    payload = { "internal_status": internal_status }
    if not test:
        (authid,authpw,server) = find_keys(key)
        
        ret = patch_obj('experiments/'+exp_id, payload, server, authid, authpw)
        print >> sys.stderr, "  * Updated encodeD %s with internal_status '%s'." % (exp_id,internal_status)
    else:
        print >> sys.stderr, "  * Would update encodeD %s with internal_status '%s'." % (exp_id,internal_status)
    return True



