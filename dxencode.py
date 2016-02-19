import os
import sys
import dxpy
import requests
import json
import re
import urlparse
import hashlib
from datetime import datetime
import time
import subprocess
import commands
import shlex

import logging

DXENCODE_VERSION = "1"

GENOME_DEFAULTS = { 'human': 'GRCh38', 'mouse': 'mm10' }
''' This the default genomes for each supported organism.'''

PRODUCTION_PROJECT = "ENCODE - Production runs"

REF_PROJECT_DEFAULT = 'ENCODE Reference Files'
''' This the default DNA Nexus project to find reference files in.'''

REF_FOLDER_DEFAULT = '/'
''' This the default folder that reference files are found in.'''


REFERENCE_FILES = {} ## Dict to cache known Reference Files
FILES = {} ## Dict to cache files
APPLETS = {} ## Dict to cache known applets

RUNS_LAUNCHED_FILE = "launchedRuns.txt"

KEYFILE = 'keypairs.json'  ## see processkey() Note this file must be in gitignore!
PRODUCTION_SERVER = 'https://www.encodeproject.org'
DEFAULT_SERVER = PRODUCTION_SERVER
S3_SERVER='s3://encode-files/'

logger = logging.getLogger('dxencode') # Callers should either use dxencode.logger or set dxencode.logger = local.logger()

def calc_md5(path):
    ''' Calculate md5 sum from file as specified by valid path name'''
    md5sum = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), ''):
            md5sum.update(chunk)
    return md5sum

SAVED_KEYS = {}

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
    return (AUTHID,AUTHPW,SERVER)
    ## TODO possibly this should return a dict
    
def encoded_post_obj(obj_type,obj_meta, SERVER, AUTHID, AUTHPW):
    ''' Posts a json object of a given type to the encoded database. '''
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
    #    print "* WARNING: request to post '%s'..." % obj_type
    #    print "Posting " + SERVER + obj_type
    #    print json.dumps(obj_meta, indent=4, sort_keys=True)
    #    print "Returned:"
    #    print json.dumps(item, indent=4, sort_keys=True)
    ##    print json.dumps(r.json(), indent=4, sort_keys=True)
    return item

def encoded_patch_obj(obj_id, obj_meta, SERVER, AUTHID, AUTHPW):
    ''' Patches a json object of a given type to the encoded database. '''
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
    #print "* request to patch %s to %s..." % (obj_id,SERVER)
    #print json.dumps(item, indent=4, sort_keys=True)
    return item
           

def encoded_post_file(filename, file_meta, SERVER, AUTHID, AUTHPW):
    ''' take a file object on local file system, post meta data and cp to AWS '''
    HEADERS = {
        'Content-type': 'application/json',
        'Accept': 'application/json',
    }

    # Handle case where file_obj is already there but file is not!
    if file_meta.get('accession') != None:
        logger.debug("Patching file object.")
        assert file_meta['status'] == "upload failed"
        file_meta['status'] = "uploading"  # Laurence will change encodeD function to take care of this
        item = encoded_patch_obj('files/'+file_meta.get('accession'),file_meta, SERVER, AUTHID, AUTHPW)
        item = encoded_post_obj('files/'+file_meta.get('accession')+"/@@upload",{}, SERVER, AUTHID, AUTHPW)
    else:
        logger.debug("Posting file object.")
        assert file_meta.get('status') == None
        item = encoded_post_obj('file',file_meta, SERVER, AUTHID, AUTHPW)

    # Could look for returned["status"] == "success"
    #returned = r.json()
    #if "status" not in returned or returned["status"] != "success" or \
    if 'upload_credentials' not in item:
        print "* ERROR: request to post %s to %s..." % (filename,SERVER)
        print json.dumps(file_meta, indent=4, sort_keys=True)
        print "* Returned..."
        print json.dumps(item, indent=4, sort_keys=True)

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
                item = encoded_patch_obj('files/'+item.get('accession'),file_meta, SERVER, AUTHID, AUTHPW)
                assert item['status'] == "upload failed"
                #raise Don't raise exception on half error... the accession needs to be added to the dx file.

    if item.get('status','uploading') != "upload failed":
        #assert item['status'] == "uploading"
        end = datetime.now()
        duration = end - start
        logger.debug("Uploaded in %.2f seconds" % duration.seconds)
            
    return item


def encoded_get(url, AUTHID=None, AUTHPW=None):
    ''' executes GET on Encoded server without without authz '''
    ##TODO possibly add try/except looking for non 4xx?
    HEADERS = {'content-type': 'application/json'}
    if AUTHID and AUTHPW:
        response = requests.get(url, auth=(AUTHID,AUTHPW), headers=HEADERS)
    else:
        response = requests.get(url, headers=HEADERS)
    return response


def enc_lookup_json(path, key='default', frame='object',must_find=False):
    '''Commonly used method to get a json object from encodeD.'''
    (AUTHID,AUTHPW,SERVER) = processkey(key)
    url = SERVER + path + '/?format=json&frame=' + frame
    #print url
    response = encoded_get(url, AUTHID, AUTHPW)
    try:
        response.raise_for_status()
        json_obj = response.json()
    except:
        if must_find:
            print "Path to json object '%s' not found." % path
            print 'Lookup failed: %s %s' % (response.status_code, response.reason)
            sys.exit(1)
        return None
    return json_obj


def env_get_current_project():
    ''' Returns the current project name for the command-line environment '''
    err, proj_name = commands.getstatusoutput('cat ~/.dnanexus_config/DX_PROJECT_CONTEXT_NAME')
    if err != 0:
        return None
    return proj_name


def project_has_folder(project, folder):
    ''' Checks for a folder in a given DX project '''
    try:
        found = project.list_folder(folder)
    except:
        return False
    return True


def folder_normalize(folder,starting=True,trailing=True):
    '''Normalizes a folder to always begin with and end with '/'.'''
    if starting:
        if not folder.startswith('/'):
            folder = '/' + folder
    else:
        if folder.startswith('/'):
            folder = folder[1:]
    if trailing:
        if not folder.endswith('/'):
            folder = folder + '/'
    else:
        if folder.endswith('/'):
            folder = folder[:-1]
    return folder

def find_folder(target_folder,project,root_folders='/',exclude_folders=["deprecated","data"]):
    '''
    Recursively attempts to find the first folder in a project and root that matches target_folder.
    The target_folder may be a nested as one/two/three but with no intervening wildcards.
    The root_folders can start with '/' but can/will grow to be nested during recursion.
    Returns full path to folder or None.
    '''
    assert len(target_folder) > 0

    # full path is easy.
    if target_folder.startswith('/') and project_has_folder(project, target_folder):
        return target_folder

    # Normalize target and root
    target_folder = folder_normalize(target_folder)
    root_folders = folder_normalize(root_folders)

    # If explicitly requesting one of the exluded folders then don't exclude it
    if exclude_folders == None:
        exclude_folders = []
    exclude_these = []
    for exclude_folder in exclude_folders:
        exclude = folder_normalize(exclude_folder)
        # Note comaprisons should always be /folder/ to /folder/ to to distinguish folder1 from redfolder10
        if root_folders.find(exclude) == -1 or target_folder.find(exclude) == -1:
            exclude_these.append(exclude)
    exclude_folders = exclude_these

    # Final normalize with target stripped of framing '/'
    target_folder = folder_normalize(target_folder,starting=False,trailing=False)

    # Because list_folder is only one level at a time, find_folder must recurse
    return rfind_folder(target_folder,project,root_folders,exclude_folders)

def rfind_folder(target_folder,project=None,root_folders='/',exclude_folders=[]):
    '''Recursive call for find_folder - DO NOT call directly.'''
    root_folders = folder_normalize(root_folders)

    for exclude_folder in exclude_folders:
        if root_folders.find(exclude_folder) != -1:
            return None
    try:
        query_folders = project.list_folder(root_folders)['folders']
    except:
        return None

    # match whole path to first target
    targets = target_folder.split('/')
    full_query = root_folders + targets[0]

    if full_query in query_folders:  # hash shortcut
        if len(targets) == 1:
            return full_query + '/' # normalized
        else:
            full_query = root_folders + target_folder # shoot for it all
            #print "- shooting [%s]" % full_query
            if project_has_folder(project, full_query):
                return full_query + '/' # normalized

    # nothing to do but recurse
    for query_folder in query_folders:
        found = rfind_folder(target_folder, project, query_folder,exclude_folders)
        if found != None:
            return found # already normalized
    return None


def find_exp_folder(project,exp_id,results_folder='/',warn=False):
    '''Returns the full path to the experiment folder if found, else None.'''
    target_folder = find_folder(exp_id,project,results_folder)
    if target_folder == None or target_folder == "":
        if warn:
            print "Unable to locate target folder (%s) for %s in project %s" % \
                                                                        (results_folder, exp_id, project.describe()['name'])
        return None
    return target_folder # already normalized


def description_from_fid(fid,properties=False):
    '''Returns file description object from fid.'''
    try:
        dxlink = FILES[fid]
    except:
        #logger.error("File %s not cached, trying id" % fid)
        dxlink = fid

    return dxpy.describe(dxlink,incl_properties=properties)


def file_handler_from_fid(fid):
    '''Returns dx file handler from fid.'''
    try:
        dxlink = FILES[fid]
    except:
        dxlink = dxpy.dxlink(fid)
    return dxpy.get_handler(dxlink)


def job_from_fid(fid):
    '''Returns job decription from fid.'''
    try:
        file_dict = description_from_fid(fid)
        job_id = file_dict["createdBy"]["job"]
        return dxpy.api.job_describe(job_id)
    except:
        return None


def file_path_from_fid(fid,projectToo=False):
    '''Returns full dx path to file from a file id.'''
    fileDict = description_from_fid(fid)
    if fileDict['folder'] == '/':
        path = '/' + fileDict['name']
    else:
        path = fileDict['folder'] + '/' + fileDict['name']
    if projectToo:
        projDict = dxpy.describe(fileDict['project'])
        path = projDict['name'] + ':' + path
    return path


def get_project(projectName, level=None):
    '''Returns the DXProject by name or errors out if not found.'''
    try:
        project = dxpy.find_one_project(name=projectName, name_mode='exact',
                                        level=level, return_handler=False)
    except:
        print "Could not find 1 and only 1 project named '"+projectName+"'."
        sys.exit(1)

    return dxpy.DXProject(project['id'])

def find_or_create_folder(project, sub_folder, root_folder='/'):
    ''' Finds or creates a sub_folder in the specified parent (root) folder'''
    folder = folder_normalize(root_folder) + folder_normalize(sub_folder,starting=False)
    if project_has_folder(project, folder):
        return folder
    else:
        logger.debug("Creating %s" % (folder))
        return project.new_folder(folder)

def get_bucket(SERVER, AUTHID, AUTHPW, f_obj):
    ''' returns aws s3 bucket and file name from encodeD file object (f_obj)'''
    #make the URL that will get redirected - get it from the file object's href property
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

    #pull out the filename
    filename = os.path.basename(o.path)

    #hack together the s3 cp url (with the s3 method instead of https)
    return filename, S3_SERVER.rstrip('/') + o.path

def copy_enc_file_to_dx(accession,proj_id,dx_folder,dx_name=None,f_obj=None,key=None):
    '''
    Finds encoded file by accession, creates s3 url and copies file to named folder,
    and adds the accession to the dx file properties, returning the dx file id or None for failure.
    If dx file name is None, the file ill have the accession based encoded name.
    MUST BE run on dx, not command-line.
    '''
    (AUTHID,AUTHPW,SERVER) = processkey(key)
    if accession != None and f_obj == None:
        url = SERVER + '/search/?type=file&accession=%s&format=json&frame=embedded&limit=all' % (accession)
        response = encoded_get(url, AUTHID, AUTHPW)
        if response == None:
            return None
        f_obj = response.json()['@graph'][0]
        if f_obj == None:
            return None
    elif accession == None and f_obj != None:
        accession = f_obj['accession']
    # If both are None or if accessions don't match then developer will hear about it!
    assert accession == f_obj['accession']

    (enc_file_name,s3_url) = get_bucket(SERVER, AUTHID, AUTHPW, f_obj)
    if enc_file_name == None or s3_url == None:
        return None
    #cp the file from the bucket
    print '> aws s3 cp %s . --quiet' %(s3_url)
    subprocess.check_call(shlex.split('aws s3 cp %s . --quiet' %(s3_url)), stderr=subprocess.STDOUT)
    subprocess.check_call(shlex.split('ls -l %s' %(enc_file_name)))
    if dx_name == None:
        dx_name = enc_file_name
        #if 'submitted_file_name' in f_obj:
        #    dx_file_name = os.path.basename(f_obj['submitted_file_name'])

    dxfile = dxpy.upload_local_file(enc_file_name,project=proj_id,folder=dx_folder, name=dx_name, \
                                    properties={ "accession": accession }, wait_on_close=True) # TODO test: wait_on_close=False
    if dxfile == None:
        return None

    return dxfile.get_id()


def move_files(fids, folder, projectId):
    '''Moves files to supplied folder.  Expected to be in the same project.'''
    for fid in fids:
        try:
            dxlink = FILES[fid]
        except:
            logger.error("File %s not in cache, trying id" % fid)
            dxlink = fid
        fileDict = dxpy.describe(dxlink) # FILES contain dxLinks
        if fileDict['project'] != projectId:
            print "ERROR: Failed to move '" + fileDict['name'] + "' as it is not in '" + \
                                                                                projectId + "'."
            sys.exit(1)
    proj = dxpy.DXProject(projectId)
    if not project_has_folder(proj, folder):
        proj.new_folder(folder,parents=True)
    proj.move(folder,fids)


def find_target_file_set(fileSet,targetFolder,project=None):
    '''Looks for a set of files in a target destination, returning the set.'''
    proj = project
    path = targetFolder
    if targetFolder.find(':') != -1:
        proj, path = targetFolder.split(':')
    if proj.find('project-') == 0:
        projId = proj
    else:
        projId = get_project(proj, level='VIEW').get_id()
    targetFids = []
    for oneFile in fileSet:
        parts = oneFile.rsplit('/',1)
        parts.reverse()
        fid = find_file(path + "/" + parts[0],projId)
        if fid != None:
            targetFids.append( fid )
    return targetFids


def find_and_copy_read_files(priors, readSet, testOnly, readToken, resultsFolder, arrayInput=False, projectId=None):
    '''Looks for read files and copies them to results folder if necessary.'''
    readFiles = []
    if len(readSet) > 0:
        readFiles = find_file_set(readSet,projectId)
        if not testOnly:
            priorReads = find_target_file_set(readSet,resultsFolder,projectId)
            if len(priorReads) == len(readFiles):
                readFiles = priorReads
            else: # any files to be moved, move all
                if readToken in priors:
                    del priors[readToken] # make sure to regenerate the combined file if necessary
                readFiles = copy_files(readFiles, projectId, resultsFolder)
        ##Incompatibility between RNA-seq and DNA-me
        if arrayInput:
            ## dna-me
            priors[readToken] = readFiles
        else:
            ## long-rna-seq back-compatibility
            if len(readFiles) > 1:        # If more than 1 file, the 'set' will need the 'concat' step.
                priors[readToken+'_set'] = readFiles
            else:                         # otherwise, the 1 file is same as 'concat' result.
                priors[readToken] = readFiles[0]

    return readFiles


def find_file_set(fileSet,projectId=None):
    '''Find all files in a set, and prints error(s) and exits if any are missing.'''
    fids = []
    if len(fileSet) > 0:
        for oneFile in fileSet:
            fid = find_file(oneFile,projectId,verbose=True)
            if fid != None:
                fids.append( fid )
        if len(fids) != len(fileSet):
            print "ERROR: expecting " + str(len(fileSet)) + " but only found " + str(len(fids)) + "."
            sys.exit(1) # verbose already gave an error message(s)

    return fids


def copy_files(fids, projectId, folder, overwrite=False):
    '''Copies array of dx file dicts to project:/folder, returning new array of dx file dicts.'''
    newFids = []
    for fid in fids:
        fileDict = dxpy.describe(FILES[fid]) # FILES contain dxLinks
        if fileDict['project'] == projectId:
            # cannot copy into the same project!!!
            # so just leave in place and pretend that we did!
            #proj = dxpy.DXProject(projectId)
            #proj.move(folder,[fid])
            newFids.append( fid )
            continue

        # Check to see if file already exists.
        alreadyThere = find_file(folder+'/'+fileDict['name'],projectId)
        if alreadyThere is None or overwrite:
            # remove what is alreadyThere?
            #if alreadyThere is not None:
            #    proj = dxpy.DXProject(projectId)
            #    proj.remove_objects([alreadyThere])
            dxFile = dxpy.get_handler(FILES[fid])
            newLink = dxpy.dxlink(dxFile.clone(projectId, folder))
        else:
            newLink = FILES(alreadyThere)
        if newLink == None:
            print "ERROR: Failed in copy of '" + fileDict['project'] + ":" + fileDict['name'] + \
                    "' to '" + projectId + ":" + folder + "'."
            sys.exit(1)
        newDict = dxpy.describe(newLink)
        FILES[newDict['id']] = newLink
        newFids.append( newDict['id'] )

    return newFids

def resolve_project(project_name, level=None):
    ''' Convert project name into DXProject object '''
    try:
        project = dxpy.find_one_project(name=project_name, name_mode='exact',
                                        level=level, return_handler=False)
    except:
        print 'Could not find 1 and only 1 project named %s; ' % format(project_name)
        sys.exit(1)

    return dxpy.DXProject(project['id'])


def get_file_link(fid, project=None):
    ''' returns a dxlink from cache or directly'''
    if not FILES.get(fid):
        FILES[fid] = dxpy.dxlink(fid,project=None)

    return FILES[fid]

def find_file(filePath,project=None,verbose=False,multiple=False, recurse=True):
    '''Using a DX style file path, find the file.'''
    proj = project
    path = filePath
    fileName = filePath
    if filePath.find(':') != -1:
        proj, path = filePath.split(':', 1)
    if path.rfind('/') != -1:
        path, fileName = path.rsplit('/', 1)
    else:
        fileName = path
        path = '/'
    if proj == None:
        if verbose:
            print "ERROR: Don't know what project to use for '" + path + "'."
        return None
    if proj.find('project-') == 0:
        projId = proj
    else:
        projId = get_project(proj, level='VIEW').get_id()
    mode = 'exact'
    if filePath.find('*') or filePath.find('?'):
        mode = 'glob'
    fileDicts = list(dxpy.find_data_objects(classname='file', folder=path, name=fileName, recurse=recurse,
                                            name_mode=mode, project=projId, return_handler=False))

    if fileDicts == None or len(fileDicts) == 0:
        #print "- Found 0 files from '" + proj + ":" + filePath + "'."
        if verbose:
            print "ERROR: Failed to find '" + proj + ":" + filePath + "'."
        return None
    elif len(fileDicts) > 1 or multiple:
        #print "- Found "+str(len(fileDict))+" files from '" + proj + ":" + filePath + "'."
        if not multiple:
            if verbose:
                print "ERROR: Found "+str(len(fileDicts))+" files when expecting 1 '" + proj + ":" + filePath + "'."
            return None
        fids = []
        for fileDict in fileDicts:
            FILES[fileDict['id']] = dxpy.dxlink(fileDict)
            fids.append( fileDict['id'] )
        return fids
    else:
        #print "- FOUND '" + proj + ":" + filePath + "'."
        FILES[fileDicts[0]['id']] = dxpy.dxlink(fileDicts[0])
        return fileDicts[0]['id']


def find_reference_file_by_name(reference_name, project_name):
    '''Looks up a reference file by name in the project that holds common tools. From Joe Dale's code.'''
    project = dxpy.find_one_project(name=project_name, name_mode='exact', return_handler=False)
    cached = '*'
    if (reference_name, project['id']) not in REFERENCE_FILES:
        found = dxpy.find_one_data_object(classname="file", name=reference_name,
                                          project=project['id'],
                                          recurse=True,
                                          zero_ok=False, more_ok=False, return_handler=True)
        REFERENCE_FILES[(reference_name, project['id'])] = found
        cached = ''

    logger.debug(cached + "Resolved %s to %s" % (reference_name, REFERENCE_FILES[(reference_name, project['id'])].get_id()))
    return dxpy.dxlink(REFERENCE_FILES[(reference_name, project['id'])])


def find_applet_by_name(applet_name, applets_project_id):
    '''Looks up an applet by name in the project that holds tools.  From Joe Dale's code.'''
    cached = '*'
    if (applet_name, applets_project_id) not in APPLETS:
        found = dxpy.find_one_data_object(classname="applet", name=applet_name,
                                          project=applets_project_id,
                                          zero_ok=False, more_ok=False, return_handler=True)
        APPLETS[(applet_name, applets_project_id)] = found
        cached = ''

    logger.debug(cached + "Resolved %s to %s" % (applet_name, APPLETS[(applet_name, applets_project_id)].get_id()))
    return APPLETS[(applet_name, applets_project_id)]

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
                print('WARNING: %s:%s Duplicate filename, ignoring.' %(exp_obj.get('accession'),file_obj.get('accession')))
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
            print("Mixed mapping for replicate %s/%s" %(experiment['accession'], rep))
            print("Paired: %s" % ([ f['accession'] for f in p ]))
            print("Unpaired: %s" % ([ f['accession'] for f in up ]))
            sys.exit(1)

    trues = len([ v for v in reps_paired.values() if v ])
    falses = len([ v for v in reps_paired.values() if not v])
    if trues and falses:
        print("Mixed mapping for replicates in experiment %s" % (experiment['accession']))
        print reps_paired
    else:
        if trues:
            return True
        elif falses:
            return False

    print "Never get here"
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
                        print "WARN: using male replacement for %s" % sex
                    sex = "male"
                organism = rep['library']['biosample']['donor']['organism']['name']
            except KeyError:
                print "Error, experiment %s replicate %s_%s missing info\n%s" % (exp_id, biorep_n, techrep_n, rep)
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

def get_exp(experiment,must_find=True,warn=False,key='default'):
    '''Returns all replicate mappings for an experiment from encoded.'''

    (AUTHID,AUTHPW,SERVER) = processkey(key)
    url = SERVER + 'experiments/%s/?format=json&frame=embedded' % experiment
    try:
        response = encoded_get(url, AUTHID, AUTHPW)
        exp = response.json()
    except:
        if must_find:
            print "Experiment %s not found." % experiment
            print "Auth %s %s may not be valid for %s." % (AUTHID, AUTHPW, url)
            #print response #.json()
            sys.exit(1)
        return None
    if exp == None or exp["status"] == "error":
        if must_find:
            print "Experiment %s not found." % experiment
            print response.json()
            sys.exit(1)
        return None

    return exp

def get_assay_type(experiment,exp=None,key='default',must_find=True,warn=False):
    '''Looks up encoded experiment's assay_type, normalized to lower case.'''
    if exp == None:
        exp = get_exp(experiment,key=key,must_find=must_find,warn=warn)
    if "assay_term_name" not in exp:
        if must_find:
            print "No 'assay_term_name' found for experiment %s." % experiment
            sys.exit(1)
        return None
    if exp["assay_term_name"] == "RNA-seq" \
    or exp["assay_term_name"] == "shRNA knockdown followed by RNA-seq" \
    or exp["assay_term_name"] == "single cell isolation followed by RNA-seq":
        #if exp["replicates"][0]["library"]["size_range"] in [">200", "300-350", "350-450"]:
        # Now more: "150-400","149-512","151-499","153-499","157-497"        
        size_range = exp["replicates"][0]["library"]["size_range"]
        if size_range.startswith('>'):
            size_range = size_range[1:]
        try:
            min_size = int(size_range.split('-')[0])
        except:
            min_size = 0
        if min_size >= 149:
            return "long-rna-seq"
        else:
            return "small-rna-seq"
    elif exp["assay_term_name"] == "whole-genome shotgun bisulfite sequencing" \
      or exp["assay_term_name"] == "shotgun bisulfite-seq assay":
        return "dna-me"
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
        print "Experiment %s has unsupported assay type of '%s'" % (experiment,exp_type)
        return None
    return exp_type

def get_full_mapping(experiment,exp=None,key='default',must_find=True,warn=False):
    '''Returns all replicate mappings for an experiment from encoded.'''

    if exp == None:
        exp = get_exp(experiment,key=key,must_find=must_find,warn=warn)

    if not exp.get('replicates') or len(exp['replicates']) < 1:
        if must_find:
            print "No replicates found in %s" % experiment
            sys.exit(1)
        return None

    return choose_mapping_for_experiment(exp,warn=warn)

def get_replicate_mapping(experiment,biorep=None,techrep=None,full_mapping=None,key='default', \
                                                                                    must_find=True):
    '''Returns replicate mappings for an experiment or specific replicate from encoded.'''
    if full_mapping == None:
        full_mapping = get_full_mapping(experiment,key=key,must_find=must_find,warn=False)

    try:
        return full_mapping[(biorep,techrep)]
    except KeyError:
        if must_find:
            print "Specified replicate: rep%s_%s could not be found in mapping of %s." % \
                ( biorep, techrep, experiment )
            print json.dumps(full_mapping,indent=4)
            sys.exit(1)
    return None

def get_reps_from_enc(exp_id, load_reads=False, exp=None, full_mapping=None, key='default'):
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
                print "Replicate has both paired(%s) and unpaired(%s) reads, quitting." % \
                    (len(mapping['paired']), len(mapping['unpaired']))
                print json.dumps(mapping,indent=4)
                sys.exit(1)                
                
            # Load read and control files, only if requested.
            run_type = None  # The files should be consistent as all single-end or paired-end, but some mapping got it wrong
            if load_reads and rep['has_reads']:
                rep['fastqs'] = { "1": [], "2": [] }
                rep['controls'] = []
                if rep['paired_end']:
                    for (p1, p2) in mapping['paired']:
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
                        rep['fastqs']['1'].append( f['accession']+".fastq.gz" )
                        if "run_type" in f:
                            if run_type == None or run_type == "single-ended":
                                run_type = f["run_type"]
                    if 'controlled_by' in mapping['unpaired']:
                        rep['controls'].append( mapping['unpaired']['controlled_by'] )
                if len(rep['controls']) == 0:
                    rep['controls'] = None
            
            # One more test because of non-standard data: single-end data reporting 'paired' == 1 !
            if rep['paired_end'] and run_type == "single-ended" and len(rep['fastqs']['2']) == 0:
                rep['paired_end'] = False
                
            reps.append( rep )

    return reps

def get_enc_file(file_acc,must_find=False,key='default'):
    '''Returns file object from encoded, given an accession.'''

    (AUTHID,AUTHPW,SERVER) = processkey(key)

    if file_acc.startswith("/files/") and file_acc.endswith('/'):
        url = SERVER + '%s?format=json&frame=embedded' % file_acc
    else:
        url = SERVER + '/files/%s/?format=json&frame=embedded' % file_acc
    try:
        response = encoded_get(url, AUTHID, AUTHPW)
        file_obj = response.json()
    except:
        if must_find:
            print "File %s not found." % file_acc
            sys.exit(1)
        return None

    return file_obj

def get_enc_exp_files(exp_obj,output_types=[],lab=None,key='default'):
    '''Returns list of file objs associated with an experiment, filtered by zero or more output_types.'''
    if not exp_obj or not exp_obj.get('files'):
        return []
    files = []
    accessions = []
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
            file_obj = get_enc_file(file_acc,key=key)
        if file_obj == None:
            continue
        #print " * Found: %s [%s] status:%s %s" % \
        #                  (file_obj['accession'],file_obj['output_type'],file_obj['status'],file_obj['submitted_file_name'])
        out_type = file_obj.get('output_type')
        if len(output_types) > 0 and (out_type == None or out_type not in output_types):
            continue
        if lab != None:
            file_lab = file_obj.get('lab')
            if file_lab != None:
                if (isinstance(file_lab,unicode) or isinstance(file_lab,str)) and file_lab != "/labs/"+lab+"/":
                    continue
                elif "name" in file_lab and file_lab["name"] != lab:
                    continue
        if file_obj.get('status') not in ["released","uploaded","uploading","in progress"]: # further restricted by caller.
            continue
        if not file_in_list(file_obj,files):
           files.extend([file_obj])
    return files

def exp_is_pe(exp,exp_files=None,rep_tech=None,server_key='default'):
    '''Determine if this experiment is expected to be 'paired-end' as opposed to 'single-end'.'''
    
    if rep_tech != None:
        reps = get_reps_from_enc(exp_id=exp.get('accession'), load_reads=False, exp=exp, full_mapping=None, key=server_key)
        for rep in reps:
            if rep.get('rep_tech') == rep_tech:
                if rep.get('paired_end') == False:
                    return False  # no more is needed!
        
    # But if rep_tech not requested or not found, OR found to be paired, then must check all fastqs
    found_fastq = False
    all_fastqs_pe = True
    
    if exp_files == None:
        exp_files = get_enc_exp_files(exp,key=server_key)
        
    for enc_f_obj in exp_files:
        if enc_f_obj.get("file_format") != "fastq":
            continue
        found_fastq = True
        if enc_f_obj.get("paired_with") == None and enc_f_obj.get("run_type") != "paired-ended":
            all_fastqs_pe = False
            
    return (found_fastq and all_fastqs_pe) 

SW_CACHE = {}
def get_sw_from_log(dxfile, regex):
    ''' given a regex and a dx file, look for the software version in the dnanexus log '''
    try:
        job_id = dxfile.describe()['createdBy']['job']
    except:
        print "Could not get job id"

    if not SW_CACHE.get(job_id+regex, {}):
        cmd = ["dx", "watch", job_id]
        log = subprocess.check_output(cmd, stderr=subprocess.STDOUT) # DEVNULL ?
        swre = re.compile(regex)
        sw = swre.findall(log)

        if not sw:
            return {}
        SW_CACHE[job_id+regex] =  {
            "software_versions":
                    [ { "software": i,
                        "version":  j }  for (i,j) in sw ]
        }
    return SW_CACHE[job_id+regex]

def create_notes(dxfile, addons={}):
    ''' creates temporary notes storage for file metadat from dxfile object '''

    description = dxfile.describe()
    notes = {
            'dx-id': description.get('id'),
            'dx-createdBy': description.get('createdBy')
    }

    notes.update(addons)
    return notes

def dx_file_get_details(fid,dxfile=None,proj_id=None):
    '''Returns dx file's details as json.'''
    if dxfile == None:
        if proj_id != None:
            dxfile = dxpy.DXFile(fid,project=proj_id)
        else:
            dxfile = file_handler_from_fid(fid)
    return dxfile.get_details()

def dx_file_get_properties(fid,dxfile=None,proj_id=None):
    '''Returns dx file's properties.'''
    if dxfile == None:
        if proj_id != None:
            dxfile = dxpy.DXFile(fid,project=proj_id)
        else:
            dxfile = file_handler_from_fid(fid)
    return dxfile.get_properties()

def dx_file_get_property(key,fid,dxfile=None,proj_id=None,return_json=False,fail_on_parse_error=True):
    '''Returns dx file's property matching 'key'.'''
    properties = dx_file_get_properties(fid,dxfile=dxfile,proj_id=proj_id)
    if not properties or key not in properties:
        return None
    if return_json:
        try:
            return json.loads(properties[key])
        except:
            try:
                return json.loads("{"+properties[key]+"}")
            except:
                print "JSON parsing failed:"
                print properties[key]
                if fail_on_parse_error:
                    sys.exit(1)
                return None
    
    return properties[key]

def dx_property_accesion_key(server):
    '''Returns the dx file propery key to use for the accession property.  Depends on the server being posted to.'''
    acc_key = "accession"
    server_key = server[server.find('/')+2:server.find('.')]# beta: "http://v25rc2.demo.encodedcc.org"
    if server_key != 'www':
        acc_key = server_key + '_accession'
    return acc_key
    
def dx_file_set_property(fid,key,value,proj_id=None,add_only=False,test=False,verbose=False):
    '''Adds/replaces key=value in a dx file's properties.
       Returns the value of the property after this operation.'''
    if proj_id != None:
        dxfile = dxpy.DXFile(fid,project=proj_id)
    else:
        dxfile = file_handler_from_fid(fid)
    properties = dxfile.get_properties()
    if verbose:
        path = '/' + dxfile.name
        if dxfile.folder != '/':
            path = dxfile.folder + path
        folder = dxfile.folder
    if key in properties:
        if properties[key] == value:
            if verbose:
                print "Note: file %s already has property '%s' set to '%s'" % (path,key,properties[key])
            return properties[key]
        elif add_only:
            if verbose:
                print "Error: file %s already has property '%s' and is not being updated from '%s'" % (path,key,properties[key])
            return properties[key]
        elif verbose:
                print "Warning: file %s has property '%s' but is being updated to '%s'" % (path,key,value)
    properties[key] = value
    if test:
        if verbose:
            print "  - Test set %s with %s='%s'" % (path,key,value)
    else:
        dxfile.set_properties(properties)
        if verbose:
            print "  - set %s with %s='%s'" % (path,key,value)
    return properties[key]
    
def umbrella_folder(folder,default,proj_name=None,exp_type=None,genome=None,annotation=None):
    '''Returns a normalized umbrella folder (that holds the experiments of a given type).'''
    if folder != default:
        return folder_normalize(folder)
            
    # No change to default, so build from parts if available
    if exp_type == None:
        return folder_normalize(folder)

    if proj_name == PRODUCTION_PROJECT:
        if exp_type == "long-rna-seq":
            folder = "/long-RNA-seq/runs/"
        elif exp_type == "small-rna-seq":
            folder = "/small-RNA-seq/runs/"
        elif exp_type == "dnase-seq":
            folder = "/DNAse-seq/runs/"
        elif exp_type == "dna-me":
            folder = "/WG Bisulfite (Methylation)/runs/"
        elif exp_type == "chip-seq":
            folder = "/ChIP-seq/runs/"
        else:
            folder = "/" + exp_type + '/runs/'
    else:
        if exp_type == "long-rna-seq":
            folder = "/lrna/"
        elif exp_type == "small-rna-seq":
            folder = "/srna/"
        elif exp_type == "dnase-seq":
            folder = "/dnase/"
        elif exp_type == "dna-me":
            folder = "/dme/"
        else:
            folder = "/" + exp_type + '/'

    if genome != None:
        folder +=  genome + '/'
        if annotation != None and genome == 'mm10':
            folder += annotation + '/'
            
    return folder


def select_alias(aliases, prefix='dnanexus:',must_find=True):
    '''returns the alias of a given prefix'''
    for alias in aliases:
        if alias.startswith(prefix):
            return alias
    if must_find:
        print "ERROR: Failed to find alias of prefix: '"+prefix+"' in:"
        print aliases
        sys.exit(1)
    
    return None


def duration_string(total_seconds,include_seconds=True):
    '''Returns a duration string.'''
    try: 
        secs = int(total_seconds)
    except:
        return None
    m, s = divmod(secs, 60)
    h, m = divmod(m, 60)
    d, h = divmod(h, 24)
    
    if include_seconds:
        if d > 0:
            return "%dd%0dh%02dm%02s" % (d, h, m, s)
        elif h > 0: 
            return    "%dh%02dm%02ds" % (h, m, s) 
        elif m > 0: 
            return        "%2dm%02ds" % (m, s)
        else: 
            return             "%2ds" % (s)
    else:
        if s >= 30:
            m +=1
        if d > 0:
            return "%dd%0dh%02dm" % (d, h, m)
        elif h > 0: 
            return     "%dh%02dm" % (h, m) 
        else: 
            return        "%2dm" % (m)


def format_duration(beg_seconds,end_seconds,include_seconds=True):
    '''Returns formatted string difference between two times in seconds.'''
    duration = end_seconds - beg_seconds
    return duration_string(duration,include_seconds)

