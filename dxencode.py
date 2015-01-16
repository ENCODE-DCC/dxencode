import os
import sys
import dxpy
import requests
import json
import re
import urlparse
import hashlib
from datetime import datetime
import subprocess

import logging

GENOME_DEFAULTS = { 'human': 'hg19', 'mouse': 'mm10' }
''' This the default genomes for each supported organism.'''

REF_PROJECT_DEFAULT = 'ENCODE Reference Files'
''' This the default DNA Nexus project to find reference files in.'''

REF_FOLDER_DEFAULT = '/'
''' This the default folder that reference files are found in.'''


REFERENCE_FILES = {} ## Dict to cache known Reference Files
FILES = {} ## Dict to cache files
APPLETS = {} ## Dict to cache known applets

RUNS_LAUNCHED_FILE = "launchedRuns.txt"

KEYFILE = 'keypairs.json'  ## see processkey() Note this file must be in gitignore!
DEFAULT_SERVER = 'https://www.encodeproject.org'
S3_SERVER='s3://encode-files/'

logger = logging.getLogger("DXENCODE")  # not sure this goes here.

def calc_md5(path):
    ''' Calculate md5 sum from file as specified by valid path name'''
    md5sum = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), ''):
            md5sum.update(chunk)
    return md5sum


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

    return (AUTHID,AUTHPW,SERVER)
    ## TODO possibly this should return a dict

def encoded_post_file(file_meta, SERVER, AUTHID, AUTHPW):
    ''' take a file object on local file system, post meta data and cp to AWS '''
    HEADERS = {
        'Content-type': 'application/json',
        'Accept': 'application/json',
    }

    r = requests.post(
        SERVER + '/file',
        auth=(AUTHID, AUTHPW),
        data=json.dumps(file_meta),
        headers=HEADERS,
    )
    try:
        r.raise_for_status()
    except:
        logger.error('Submission failed: %s %s' % (r.status_code, r.reason))
        logger.error(r.text)
        raise
    item = r.json()['@graph'][0]

    ####################
    # POST file to S3

    creds = item['upload_credentials']
    env = os.environ.copy()
    env.update({
        'AWS_ACCESS_KEY_ID': creds['access_key'],
        'AWS_SECRET_ACCESS_KEY': creds['secret_key'],
        'AWS_SECURITY_TOKEN': creds['session_token'],
    })

    logger.debug("Uploading file.")
    start = datetime.now()
    try:
        subprocess.check_call(['aws', 's3', 'cp', file_meta['submitted_file_name'], creds['upload_url']], env=env)
        end = datetime.now()
        duration = end - start
        logger.debug("Uploaded in %.2f seconds" % duration.seconds)
    except:
        logger.debug("Upload failed")

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


def project_has_folder(project, folder):
    ''' Checks for a folder in a given DX project '''
    try:
        found = project.list_folder(folder)
    except:
        return False
    return True


def find_folder(target_folder,project,root_folders='/',exclude_deprecated=True):
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
    if target_folder.endswith('/'):
        target_folder = target_folder[:-1]
    if root_folders[0] != '/':
       root_folders = '/' + root_folders
       
    # If explicitly requesting deprecated, don't exclude it
    if exclude_deprecated:
        if root_folders.find('/deprecated') != -1 or target_folder.find('/deprecated') != -1:
            exclude_deprecated=False

    # Because list_folder is only one level at a time, find_folder must recurse
    return rfind_folder(target_folder,project,root_folders,exclude_deprecated)  

def rfind_folder(target_folder,project=None,root_folders='/',exclude_deprecated=True):
    '''Recursive call for find_folder - DO NOT call directly.'''
    if exclude_deprecated and root_folders.find('/deprecated') != -1:
        return None
    try:
        query_folders = project.list_folder(root_folders)['folders']
    except:
        return None
    
    # Normalize    
    if root_folders.endswith('/'):
        root_folders = root_folders[-1]
        
    # match whole path to first target
    targets = target_folder.split('/')
    full_query = root_folders + '/' + targets[0] 
    #print "Full query [%s]" % full_query

    if full_query in query_folders:  # hash shortcut
        if len(targets) == 1:
            return full_query
        else:
            full_query = root_folders + '/' + target_folder # shoot for it all
            #print "- shooting [%s]" % full_query
            if project_has_folder(project, full_query):
                return full_query

    # nothing to do but recurse
    for query_folder in query_folders:
        found = rfind_folder(target_folder, project, query_folder,exclude_deprecated)
        if found != None:
            return found
    return None      

def file_path_from_fid(fid,projectToo=False):
    '''Returns full dx path to file from a file id.'''
    try:
        dxlink = FILES[fid]
    except:
        logger.error("File %s not cached, trying id" % fid)
        dxlink = fid

    fileDict = dxpy.describe(dxlink) # FILES contain dxLinks
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
    folder = root_folder+sub_folder
    logger.debug("Creating %s (%s)" % (folder, root_folder))
    if project_has_folder(project, folder):
        return folder
    else:
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
            targetFids += [ fid ]
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
                fids += [ fid ]
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
            newFids += [ fid ]
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
        newFids += [ newDict['id'] ]

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
            fids += [ fileDict['id'] ]
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

def filenames_in(files=None):
    if not len(files):
        return []
    else:
        return [f.get('submitted_file_name') for f in files]

def files_to_map(exp_obj):
    if not exp_obj or not exp_obj.get('files'):
        return []
    else:
        files = []
        for file_obj in exp_obj.get('files'):
            if (file_obj.get('output_type') == 'reads' or file_obj.get('output_type') == 'raw data') and \
               file_obj.get('file_format') == 'fastq' and \
               file_obj.get('submitted_file_name') not in filenames_in(files):
               files.extend([file_obj])
            elif file_obj.get('submitted_file_name') in filenames_in(files):
                logger.warning('%s:%s Duplicate filename, ignoring.' %(exp_obj.get('accession'),file_obj.get('accession')))
                return []
        return files

def replicates_to_map(experiment, files):
    if not files:
        return []
    else:
        reps_with_files = set([ f['replicate']['uuid'] for f in files if f.get('replicate') ])
        return [ r for r in experiment['replicates'] if r['uuid'] in reps_with_files ]

def choose_mapping_for_experiment(experiment,warn=True):
    ''' for a given experiment object, fully embedded, return experimental info needed for mapping
        returns an dict keyed by [biological_rep][technical_rep]
        with information for mapping (sex, organism, paired/unpaired files, library id)
    '''

    exp_id = experiment['accession']
    files = files_to_map(experiment)
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
                        logging.warning('%s:%s could not find mate' %(experiment.get('accession'), file_object.get('accession')))
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
                logging.warning('%s: leftover file(s) %s' % (exp_id, rep_files))
    elif warn:
        logging.warning('%s: No files to map' % exp_id)
    return mapping

def get_exp(experiment,must_find=True,warn=False):
    '''Returns all replicate mappings for an experiment from encoded.'''
    
    (AUTHID,AUTHPW,SERVER) = processkey('default')
    url = SERVER + 'experiments/%s/?format=json&frame=embedded' % experiment
    response = encoded_get(url, AUTHID, AUTHPW)
    exp = response.json()

    if not exp.get('replicates') or len(exp['replicates']) < 1:
        if must_find:
            print "No replicates found in %s\n%s" % ( experiment, exp )
            sys.exit(1)
        return None

    return exp

def get_full_mapping(experiment,exp=None,must_find=True,warn=False):
    '''Returns all replicate mappings for an experiment from encoded.'''
    
    if exp == None:
        exp = get_exp(experiment,must_find=must_find,warn=warn)

    return choose_mapping_for_experiment(exp,warn=warn)

def get_replicate_mapping(experiment,biorep=None,techrep=None,full_mapping=None,must_find=True):
    '''Returns replicate mappings for an experiment or specific replicate from encoded.'''
    if full_mapping == None:
        full_mapping = get_full_mapping(experiment,must_find=True,warn=False)
        
    try:
        return full_mapping[(biorep,techrep)]
    except KeyError:
        if must_find:
            print "Specified replicate: rep%s_%s could not be found in mapping of %s." % \
                ( biorep, techrep, experiment )
            print json.dumps(full_mapping,indent=4)
            sys.exit(1)
    return None

def load_fastqs_from_mapping(load_to, mapping, controls=False):
    '''
    Resolves fastq file names from mapping and puts them into load_to dict (hint psv).
    Note: load_to['fastqs'] will always be { "1":[],"2":[]} though '2' is an empty set for unpaired
    Returns True or False for paired end.  Failures exit.
    '''

    # Paired ends?  Read files?
    if mapping['unpaired'] and not mapping['paired']:
        paired_end = False
    elif mapping['paired'] and not mapping['unpaired']:
        paired_end = True
    elif not mapping['unpaired'] and not mapping['paired']:
        print "Replicate has no reads either paired or unpaired"
        print json.dumps(mapping,indent=4)
        sys.exit(1)
    else:
        print "Replicate has both paired(%s) and unpaired(%s) reads, quitting." % \
            (len(mapping['paired'], len(mapping['unpaired'])))
        print json.dumps(mapping,indent=4)
        sys.exit(1)
    if controls:
        load_to['controls'] = []
    load_to['fastqs'] = { "1": [], "2": [] }
    if paired_end:
        for (p1, p2) in mapping['paired']:
            load_to['fastqs'][p1['paired_end']].append(p1['accession']+".fastq.gz")
            if p2 != None and 'paired_end' in p2:
                load_to['fastqs'][p2['paired_end']].append(p2['accession']+".fastq.gz")
            if controls:
                if 'controlled_by' in p1:
                    load_to['controls'] += p1['controlled_by']
                if p2 != None and 'controlled_by' in p2:
                    load_to['controls'] += p2['controlled_by']
    else:
        load_to['fastqs']['1'] = [ f['accession']+".fastq.gz" for f in mapping['unpaired'] ]
        if controls:
            if 'controlled_by' in mapping['unpaired']:
                load_to['controls'] += mapping['unpaired']['controlled_by']

    return paired_end

def common_variables(args,results_folder_default,fastqs=True,controls=False):
    '''Initializes dict with common variables from args and encoded.'''
    
    cv = {}
    cv['project']    = args.project
    cv['experiment'] = args.experiment

    # expecting either commbined-replicates or biological and technical replicate
    if 'cr' in args and args.cr != None:
        if len(args.cr) != 2:
            print "Specify which two replicates to compare (e.g. '1_2 2')."
            sys.exit(1)
        cv['combined'] = True
        cv['reps'] = {}
        for rep_tech in args.cr:
            # Normalize rep_tech string
            if '_' not in rep_tech:
                rep_tech = rep_tech + '_1' # default to tech_rep 1
            if not rep_tech.startswith('rep'):
                rep_tech = 'rep' + rep_tech
            # parse out biological and technical replicates
            br_tr = rep_tech[3:]
            (br,tr) = br_tr.split('_')
            if 'a' in cv['reps']:
                cv['reps']['b'] = { 'br': int(br), 'tr': int(tr) }
                cv['reps']['b']['rep_tech'] = rep_tech
            else:
                cv['reps']['a'] = { 'br': int(br), 'tr': int(tr) }
                cv['reps']['a']['rep_tech'] = rep_tech
                
        if cv['reps']['a']['rep_tech'] == cv['reps']['b']['rep_tech']:
            print "Specify different replicates to compare (e.g. 'rep1_1 rep2_1')."
            sys.exit(1)
        args.br = cv['reps']['a']['br']
        args.tr = cv['reps']['a']['tr']
        cv['rep_tech'] = 'reps' + cv['reps']['a']['rep_tech'][3:] + \
                            '-' + cv['reps']['b']['rep_tech'][3:]
    else:
        cv['combined'] = False
        if args.br == 0:
            print "Must specify either --biological-replicate or --compare-replicates." 
            sys.exit(1)
        cv['reps'] = {'a': { 'br': args.br, 'tr': args.tr} }
        cv['reps']['a']['rep_tech'] = 'rep' + str(args.br) + '_' + str(args.tr)

    full_mapping = get_full_mapping(cv['experiment'])
    for cv_rep in cv['reps'].keys():
        rep = cv['reps'][cv_rep]
        # TODO get rep_mapping once and then subset "mapping" multiple times
        mapping = get_replicate_mapping(cv['experiment'],rep['br'],rep['tr'], full_mapping)

        # Not a common var but convenient and cheap
        rep['library_id'] = mapping['library']
        rep['replicate_id'] = mapping['replicate_id']
        
        # TODO add enough info that individual pipelines can varify the experiment matches pipeline
        
        if fastqs:
            rep['paired_end'] = load_fastqs_from_mapping(rep,mapping,controls)
            # Non-file app inputs
            rep['rootR1'] = cv['experiment'] + rep['rep_tech'] + '_concatR1'
            if rep['paired_end']: 
                rep['rootR2'] = cv['experiment'] + rep['rep_tech'] + '_concatR2'

        if cv_rep == 'a':
            # Only supported genomes
            if mapping['organism'] in GENOME_DEFAULTS:
                cv['genome'] = GENOME_DEFAULTS[mapping['organism']]
            else:
                print "Organism %s not currently supported" % mapping['organism']
                sys.exit(1)

            cv['gender'] = mapping['sex']

    # Paired ends?
    if cv['combined'] and cv['reps']['a']['paired_end'] != cv['reps']['a']['paired_end']:
        print "Replicates are expected to be both paired-end or single-end!  Check encoded."
        sys.exit(1)
    cv['paired_end'] = cv['reps']['a']['paired_end']

    # Default locations (with adjustments)
    cv['refLoc'] = args.refLoc
    if cv['refLoc'] == REF_FOLDER_DEFAULT:
        cv['refLoc'] = REF_FOLDER_DEFAULT + cv['genome'] + '/'
    if not cv['refLoc'].endswith('/'):
        cv['refLoc'] += '/' 
    cv['resultsLoc'] = args.resultsLoc
    if cv['resultsLoc'] == results_folder_default:
        if cv['genome'] == 'mm10':
            cv['resultsLoc'] = results_folder_default + cv['genome'] + '/' + cv['annotation'] + '/'
        else:
            cv['resultsLoc'] = results_folder_default + cv['genome'] + '/'
    if not cv['resultsLoc'].endswith('/'):
        cv['resultsLoc'] += '/' 
    cv['resultsFolder'] = cv['resultsLoc'] + cv['experiment'] + '/'
    cv['reps']['a']['resultsFolder'] = cv['resultsLoc'] + cv['experiment'] + '/' + \
                                                          cv['reps']['a']['rep_tech'] + '/'
    if cv['combined']:
        cv['reps']['b']['resultsFolder'] = cv['resultsLoc'] + cv['experiment'] + '/' + \
                                                              cv['reps']['b']['rep_tech'] + '/'

    return cv

def find_prior_results(pipe_path,steps,results_folder,file_globs,proj_id):
    '''Looks for all result files in the results folder.'''
    priors = {}
    for step in pipe_path:
        for fileToken in steps[step]['results'].keys():
            fid = find_file(results_folder + file_globs[fileToken],proj_id,recurse=False)
            if fid != None:
                priors[fileToken] = fid
    return priors

def finding_rep_inputs_and_priors(psv,steps,globs,project,test):
    '''Finds the inputs and priors for a run.'''
    
    print "Checking for prior results..."
    proj_id = project.get_id()

    # NOTE: priors is a dictionary of fileIds that will be used to determine stepsToDo
    #       and fill in inputs to workflow steps
    for rep in psv['reps'].values():
        if not test:
            if not project_has_folder(project, rep['resultsFolder']):
                project.new_folder(rep['resultsFolder'],parents=True)
        rep['priors'] = find_prior_results(rep['path'],steps,rep['resultsFolder'],globs,proj_id)

    print "Checking for input files..."
    # Find all reads files and move into place
    # TODO: files could be in: dx (usual), remote (url e.g.https://www.encodeproject.org/...
    #       or possibly local, Currently only DX locations are supported.
    for rep in psv['reps'].values():
        rep['inputs'] = {}
        reads_token = 'reads'
        if psv['paired_end']:
            reads_token += '1'
        rep['inputs']['Reads1'] = find_and_copy_read_files(rep['priors'], rep['fastqs']['1'], \
                                           test, reads_token, rep['resultsFolder'], False, proj_id)
        # Note: rep['fastqs']['2'] and rep['inputs']['Reads2'] will be empty on single-end
        rep['inputs']['Reads2'] = find_and_copy_read_files(rep['priors'], rep['fastqs']['2'], \
                                            test, 'reads2', rep['resultsFolder'], False, proj_id)

def find_all_ref_files(psv,find_references_function):
    '''Locates all reference files based upon organism and gender.'''
    
    print "Looking for reference files..."
    for rep in psv['reps'].values():
        # Need multiple copies of control files because they are put into priors!
        find_references_function(rep['priors'],psv) # pipeline specific
    if psv['combined']:
        find_references_function(psv['priors'],psv)


def determine_steps_to_run(pipe_path, steps, priors, deprecate, proj_id, force=False, verbose=False):
    '''Determine what steps need to be done, base upon prior results.'''
    will_create = []
    steps_to_run = []
    for step in pipe_path:
        # Force will include the first step with all its inputs
        # This should avoid forcing concat if it isn't needed
        #
        if force:
            inputs = steps[step]['inputs'].keys()
            count = 0
            for input in inputs:
                if input in priors:
                    count += 1
            if count == len(inputs):
                steps_to_run += [ step ]
                if verbose:
                    print "- Adding step '"+step+"' because of force flag."
        if step not in steps_to_run:
            results = steps[step]['results'].keys()
            for result in results:
                if result not in priors:
                    steps_to_run += [ step ]
                    if verbose:
                        print "- Adding step '"+step+"' because prior '"+result+"' was not found."
                    break
        # If results are there but inputs are being recreated, then step must be rerun
        if step not in steps_to_run:
            inputs = steps[step]['inputs'].keys()
            for inp in inputs:
                if inp in will_create:
                    steps_to_run += [ step ]
                    if verbose:
                        print "- Adding step '"+step+"' due to prior step dependency."
                    break
        # Any step that is rerun, will cause prior results to be deprecated
        # NOTE: It is necessary to remove from 'priors' so succeeding steps are rerun
        # NOTE: It is also important to move prior results out of target folder to avoid confusion!
        if step in steps_to_run:
            results = steps[step]['results'].keys()
            for result in results:
                will_create += [ result ]
                if result in priors:
                    deprecate += [ priors[result] ]
                    del priors[result]
                    # if results are in folder, then duplicate files cause a problem!
                    # So add to 'deprecate' to move or remove before launching

    # Now make sure the steps can be found, and error out if not.
    for step in steps_to_run:
        app = steps[step]['app']
        dxApp = dxpy.find_data_objects(classname='file', name=app, name_mode='exact',
                                                         project=proj_id, return_handler=False)
        if dxApp == None:
            print "ERROR: failure to locate app '"+app+"'!"
            sys.exit(1)

    return steps_to_run


def determine_steps_needed(psv, replicate_steps, combined_steps, proj_id, force=False):
    '''Determine steps needed for replicate(s) and combined, base upon prior results.'''
    
    print "Determining steps to run..."
    # NOTE: stepsToDo is an ordered list of steps that need to be run
    for rep in psv['reps'].values():
        rep['deprecate'] = [] # old results will need to be moved/removed if step is rerun
        rep['stepsToDo'] = determine_steps_to_run(rep['path'], replicate_steps, \
                                     rep['priors'], rep['deprecate'], proj_id, force=force)
    if psv['combined']:
        psv['deprecate'] = [] # old results will need to be moved/removed if step is rerun
        psv['stepsToDo'] = determine_steps_to_run(psv['path'], combined_steps, \
                                     psv['priors'], psv['deprecate'], proj_id, force=force)


def check_run_log(results_folder,proj_id,verbose=False):
    '''Checks for currently running jobs and will exit if found.'''
    run_log_path = results_folder + '/' + RUNS_LAUNCHED_FILE
    log_fids = find_file(run_log_path,proj_id,multiple=True,recurse=False)
    if log_fids == None:
        if verbose:
            print "  No prior jobs launched."
    else:
        # NOTE: Appending to the one file, but just in case handle multiple files.
        for fid in log_fids:
            with dxpy.open_dxfile(fid) as fd:
                for line in fd:
                    #print "Looking for job ["+line+"]"
                    run_id = line.split(None,1)
                    if not run_id[0].startswith('analysis-'):
                        continue
                    analysis = dxpy.DXAnalysis(dxid=run_id[0])
                    if analysis == None:
                        continue
                    state = analysis.describe()['state']
                    # states I have seen: in_progress, terminated, done, failed
                    if state not in [ "done", "failed", "terminated" ]:
                        msg="Exiting: Can't launch because prior run ["+run_id[0]+"] "
                        if len(run_id) > 1:
                            msg+="("+run_id[1]+") "
                        msg+= "has not finished (currently '"+state+"')."
                        print msg
                        sys.exit(1)
                    elif verbose:
                        msg="  Prior run ["+run_id[0]+"] "
                        if len(run_id) > 1:
                            msg+="("+run_id[1]+") "
                        msg+= "is '"+state+"'."
                        print msg

def log_this_run(run_id,results_folder,proj_id):
    '''Adds a runId to the runsLaunched file in resultsFolder.'''
    # NOTE: DX manual lies?!  Append not possible?!  Then write new/delete old
    run_log_path = results_folder + '/' + RUNS_LAUNCHED_FILE
    old_fids = find_file(run_log_path,proj_id,multiple=True,recurse=False)
    new_fh = dxpy.new_dxfile('a',project=proj_id,folder=results_folder,name=RUNS_LAUNCHED_FILE)
    new_fh.write(run_id+' started:'+str(datetime.now())+'\n')
    if old_fids is not None:
        for old_fid in old_fids:
            with dxpy.open_dxfile(old_fid) as old_fh:
                for old_run_id in old_fh:
                    new_fh.write(old_run_id+'\n')
        proj = dxpy.DXProject(proj_id)
        proj.remove_objects(old_fids)
    new_fh.close()

def create_workflow(psv, run, proj_id, app_proj_id=None,test=False,template=False):
    '''
    This function will populate a workflow for the steps in run['stepsToDo'] and return 
    the worklow unlaunched.   It relies on steps dict which contains input and output requirements,
    pvs (pipeline specific variables) dictionary and run (run specific dictionary),
    which contains input and previous results already in results dir
    '''

    if len(run['stepsToDo']) < 1:
        return None
    if app_proj_id == None:
        app_proj_id = proj_id

    # create a workflow object
    if not test:
        wf = dxpy.new_dxworkflow(title=run['name'],name=run['name'],folder=run['resultsFolder'],
                                                project=proj_id,description=psv['description'])

    # NOTE: prevStepResults dict contains links to result files to be generated by previous steps
    prevStepResults = {}
    steps = run['steps']
    for step in run['stepsToDo']:
        appName = steps[step]['app']
        app = find_applet_by_name(appName, app_proj_id)
        appInputs = {}
        # file inputs
        for fileToken in steps[step]['inputs'].keys():
            appInp = steps[step]['inputs'][fileToken]
            if fileToken in prevStepResults:
                appInputs[ appInp ] = prevStepResults[fileToken]
            elif fileToken in run['priors']:
                if isinstance(run['priors'][fileToken], list):
                    appInputs[ appInp ] = []
                    for fid in run['priors'][fileToken]:
                        appInputs[ appInp ] += [ FILES[fid] ]
                else:
                    appInputs[ appInp ] = FILES[ run['priors'][fileToken] ]
            elif not template:
                print "ERROR: step '"+step+"' can't find input '"+fileToken+"'!"
                sys.exit(1)
        # Non-file app inputs
        if 'params' in steps[step]:
            for param in steps[step]['params'].keys():
                appParam = steps[step]['params'][param]
                if param in run:
                    appInputs[ appParam ] = run[param]
                elif param in psv:
                    appInputs[ appParam ] = psv[param]
                elif not template:
                    print "ERROR: step '"+step+"' unable to locate '"+param+ \
                                                        "' in pipeline specific variables (psv)."
                    sys.exit(1)
        # Add wf stage
        if not test:
            stageId = wf.add_stage(app, stage_input=appInputs, folder=run['resultsFolder'])
        # outputs, which we will need to link to
        for fileToken in steps[step]['results'].keys():
            appOut = steps[step]['results'][fileToken]
            if test:
                prevStepResults[ fileToken ] = 'fake-for-testing'
            else:
                prevStepResults[ fileToken ] = dxpy.dxlink({ 'stage': stageId,'outputField': appOut })

    if test:
        return None
    else:
        return wf


def build_a_step(applet, file_globs, proj_id):
    ''' create input object for a step and extends the file_globs dict as appropriate.'''
    try:
        dxapp = dxpy.find_one_data_object(classname="applet", name=applet, project=proj_id,
                                          zero_ok=False, more_ok=False,describe=True)
    except IOError:
        print "Cannot find applet '"+applet+"' in "+proj_id
        sys.exit(1)

    params = {}
    inputs = {}
    results = {}
    inps = dxapp['describe'].get('inputSpec') or []
    for inp in inps:
        in_name = inp['name'].encode('ascii','ignore')
        if inp['class'] == 'file' or inp['class'] == 'array:file':
            inputs[in_name] = in_name
        else:
            params[in_name] = in_name

    outs = dxapp['describe'].get('outputSpec') or []
    for out in outs:
        if out['class'] == 'file' or out['class'] == 'array:file':
            out_name = out['name'].encode('ascii','ignore')
            results[out_name] = out_name
            if out_name not in file_globs and 'patterns' in out:
                file_globs[out_name] = '/'+out['patterns'][0].encode('ascii','ignore')
        else:
            pass
            # TODO not sure what to do with these

    return {
        'app': applet,
        'params': params,
        'inputs': inputs,
        'results': results
    }

def build_simple_steps(pipe_path, proj_id, verbose=False):
    '''
    Builds dict of steps for the apps in the pipeline and a dict of file_globs for look up.
    Only works for pipelines where every step in pipe_path is a distinct app,
    and each result glob is uniq for all pipeline results.
    '''
    steps = {}
    file_globs = {}
    for step in pipe_path:
        steps[step] = build_a_step(step, file_globs, proj_id)
    if verbose:
        print "STEPS = "
        print json.dumps(steps,indent=4)
        print "FILE_GLOBS = "
        print json.dumps(file_globs,indent=4)

    return [ steps, file_globs ]


def report_run_plans(psv, run):
    '''Report the plans before executing them.'''

    print "Running: "+run['title']
    if 'subTitle' in run:
        print "         "+run['subTitle']
    for input_type in sorted( run['inputs'].keys() ):
        if len(run['inputs'][input_type]) > 0: 
            print "- " + input_type + ":"
            for fid in run['inputs'][input_type]:
                print "  " + file_path_from_fid(fid)
    print "- Reference files:"
    for token in psv['ref_files']:
        print "  " + file_path_from_fid(run['priors'][token],True)
    print "- Results written to: " + psv['project'] + ":" +run['resultsFolder']
    if len(run['stepsToDo']) == 0:
        print "* All expected results are in the results folder, so there is nothing to do."
        print "  If this experiment/replicate needs to be rerun, then use the --force flag to "
        print "  rerun all steps; or remove suspect results from the folder before relaunching."
        sys.exit(0)
    else:
        print "- Steps to run:"
        for step in run['path']:
            if step in run['stepsToDo']:
                print "  * "+run['steps'][step]['app']+" will be run"
            else:
                if not step.find('concat') == 0:
                    print "    "+run['steps'][step]['app']+" has already been run"

    if len(run['deprecate']) > 0:
        deprecated = run['resultsFolder']+"/deprecated/"
        if run['resultsFolder'].endswith('/'):
            deprecated = run['resultsFolder']+"deprecated/"
        print "Will move "+str(len(run['deprecate']))+" prior result file(s) to '" + deprecated+"'."
        for fid in run['deprecate']:
            print "  " + file_path_from_fid(fid)


def launchPad(wf,proj_id,psv,run=False):
    '''Launches or just advertises preassembled workflow.'''
    if wf == None:
        print "ERROR: failure to assemble workflow!"
        sys.exit(1)

    if run:
        print "Launch sequence initiating..."
        wf_run = wf.run({}, project=proj_id)
        if wf_run == None:
            print "ERROR: failure to lift off!"
            sys.exit(1)
        else:
            print "  We have liftoff!"
            wf_dict = wf_run.describe()
            log_this_run(wf_dict['id'],psv['resultsFolder'],proj_id)
            print "  Launched " + wf_dict['id']+" as '"+wf.name+"'"
    else:
        print "Workflow '" + wf.name + "' has been assembled in "+psv['resultsFolder'] + \
                                                                    ". Manual launch required."


def report_build_launch(psv, run, proj_id, test=True, launch=False):
    '''
    Handles everything from reporting to building workflow to running, using standardized 
    methods.  This function should finish a launcher using the psv and run dict already created.
    '''
    # Report the plans
    report_run_plans(psv, run)

    print "Checking for currently running analyses..."
    check_run_log(run['resultsFolder'], proj_id, verbose=True)

    if len(run['deprecate']) > 0 and not test:
        deprecated = run['resultsFolder']+"deprecated/"
        print "Moving "+str(len(run['deprecate']))+" prior result file(s) to '"+deprecated+"'..."
        move_files(run['deprecate'],deprecated,proj_id)

    if test:
        print "Testing workflow assembly..."
    else:
        print "Assembling workflow..."
    wf = create_workflow(psv, run, proj_id, test=test,template=(not launch))

    # Exit if test only
    if test:
        print "TEST ONLY - exiting."
        sys.exit(0)

    # Roll out to pad and possibly launch
    launchPad(wf,proj_id,run,launch)


SW_CACHE = {}
def get_sw_from_log(dxfile, regex):
    ''' given a regex and a dx file, look for the software version in the dnanexus log '''
    try:
        job_id = dxfile.describe()['createdBy']['job']
    except:
        print "Could not get job id"

    if not SW_CACHE.get(job_id, {}):
        cmd = ["dx", "watch", job_id]
        log = subprocess.check_output(cmd)
        swre = re.compile(regex)
        sw = swre.findall(log)

        if not sw:
            return {}
        SW_CACHE[job_id] =  {
            "software_versions":
                    [ { "software": i,
                        "version":  j }  for (i,j) in sw ]
        }
    return SW_CACHE[job_id]

def create_notes(dxfile, addons={}):
    ''' creates temporary notes storage for file metadat from dxfile object '''

    description = dxfile.describe()
    notes = {
            'dx-id': description.get('id'),
            'dx-createdBy': description.get('createdBy')
    }

    notes.update(addons)
    return notes


