import os, sys, json
import subprocess, commands
import hashlib, re
import dxpy
#import shlex

import logging

DX_VERSION = "1"

GENOME_DEFAULTS = { 'human': 'GRCh38', 'mouse': 'mm10' }
''' This the default genomes for each supported organism.'''

PRODUCTION_PROJECT = "ENCODE - Production runs"

REF_PROJECT_DEFAULT = 'ENCODE Reference Files'
''' This the default DNA Nexus project to find reference files in.'''

REF_FOLDER_DEFAULT = '/'
''' This the default folder that reference files are found in.'''

INTERNAL_STATUS_BLOCKS = ["requires lab review", "unrunnable"]
'''Experients with these internal_statuses should not be assembled or launched.'''

REFERENCE_FILES = {} ## Dict to cache known Reference Files
FILES = {} ## Dict to cache files
APPLETS = {} ## Dict to cache known applets

RUNS_LAUNCHED_FILE = "launchedRuns.txt"
    
def clear_cache():
    '''Empties all cache'''
    global REFERENCE_FILES
    global FILES
    global APPLETS
    REFERENCE_FILES = {} ## Dict to cache known Reference Files
    FILES = {} ## Dict to cache files
    APPLETS = {} ## Dict to cache known applets

def calc_md5(path):
    ''' Calculate md5 sum from file as specified by valid path name'''
    md5sum = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), ''):
            md5sum.update(chunk)
    return md5sum

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
        query_folders = project.list_folder(root_folders,only='folders')['folders']
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
            #print >> sys.stderr, "- shooting [%s]" % full_query
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
            print >> sys.stderr, "Unable to locate target folder (%s) for %s in project %s" % \
                                                                        (results_folder, exp_id, project.describe()['name'])
        return None
    return target_folder # already normalized


def find_replicate_folders(project,exp_folder,verbose=False):
    '''Returns a sorted list of replicate folders which are sub-folders of an exp folder.'''
    # normalize
    try:
        sub_folders = project.list_folder(exp_folder,only='folders')['folders']
    except:
        if verbose:
            print >> sys.stderr, "No subfolders found for %s" % exp_folder
        return []

    rep_folders = []
    for path in sorted(sub_folders):
        folder = path.split('/')[-1]
        if folder.startswith('rep'):
            if len(folder[3:].split('_')) == 2:
                rep_folders.append( folder )
    if verbose:
        print >> sys.stderr, "Replicate folders:"
        print >> sys.stderr, json.dumps(rep_folders,indent=4)
    return rep_folders


def description_from_fid(fid,properties=False):
    '''Returns file description object from fid.'''
    try:
        dxlink = FILES[fid]
    except:
        #print >> sys.stderr, "File %s not cached, trying id" % fid)
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
        print >> sys.stderr, "Creating %s" % (folder)
        return project.new_folder(folder)

def move_files(fids, folder, projectId):
    '''Moves files to supplied folder.  Expected to be in the same project.'''
    for fid in fids:
        try:
            dxlink = FILES[fid]
        except:
            #print >> sys.stderr, "File %s not in cache, trying id" % fid
            dxlink = fid
        fileDict = dxpy.describe(dxlink) # FILES contain dxLinks
        if fileDict['project'] != projectId:
            print >> sys.stderr, "ERROR: Failed to move '" + fileDict['name'] + "' as it is not in '" + \
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
        else:
            if verbose:
                print " Found "+str(len(fileDicts))+" files for '" + proj + ":" + filePath + "'."
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
    cached = '* '
    if (reference_name, project['id']) not in REFERENCE_FILES:
        found = dxpy.find_one_data_object(classname="file", name=reference_name,
                                          project=project['id'],
                                          recurse=True,
                                          zero_ok=False, more_ok=False, return_handler=True)
        REFERENCE_FILES[(reference_name, project['id'])] = found
        cached = ''

    #print >> sys.stderr, cached + "Resolved %s to %s" % \
    #                                                (reference_name, REFERENCE_FILES[(reference_name, project['id'])].get_id())
    return dxpy.dxlink(REFERENCE_FILES[(reference_name, project['id'])])


def find_applet_by_name(applet_name, applets_project_id):
    '''Looks up an applet by name in the project that holds tools.  From Joe Dale's code.'''
    cached = '* '
    if (applet_name, applets_project_id) not in APPLETS:
        found = dxpy.find_one_data_object(classname="applet", name=applet_name,
                                          project=applets_project_id,
                                          zero_ok=False, more_ok=False, return_handler=True)
        APPLETS[(applet_name, applets_project_id)] = found
        cached = ''

    #print >> sys.stderr, cached + "Resolved %s to %s" % (applet_name, APPLETS[(applet_name, applets_project_id)].get_id())
    return APPLETS[(applet_name, applets_project_id)]

SW_CACHE = {}
def get_sw_from_log(dxfile, regex):
    ''' given a regex and a dx file, look for the software version in the dnanexus log '''
    try:
        job_id = dxfile.describe()['createdBy']['job']
    except:
        print >> sys.stderr, "Could not get job id"

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

def file_get_details(fid,dxfile=None,proj_id=None):
    '''Returns dx file's details as json.'''
    if dxfile == None:
        if proj_id != None:
            dxfile = dxpy.DXFile(fid,project=proj_id)
        else:
            dxfile = file_handler_from_fid(fid)
    return dxfile.get_details()

def file_get_properties(fid,dxfile=None,proj_id=None):
    '''Returns dx file's properties.'''
    if dxfile == None:
        if proj_id != None:
            dxfile = dxpy.DXFile(fid,project=proj_id)
        else:
            dxfile = file_handler_from_fid(fid)
    return dxfile.get_properties()

def file_get_property(key,fid,dxfile=None,proj_id=None,return_json=False,fail_on_parse_error=True):
    '''Returns dx file's property matching 'key'.'''
    properties = file_get_properties(fid,dxfile=dxfile,proj_id=proj_id)
    if not properties or key not in properties:
        return None
    if return_json:
        try:
            return json.loads(properties[key])
        except:
            try:
                return json.loads("{"+properties[key]+"}")
            except:
                print >> sys.stderr, "JSON parsing failed:"
                print >> sys.stderr, properties[key]
                if fail_on_parse_error:
                    sys.exit(1)
                return None
    
    return properties[key]

def property_accesion_key(server):
    '''Returns the dx file propery key to use for the accession property.  Depends on the server being posted to.'''
    acc_key = "accession"
    if server.find('/') == -1:
        server_key = server
    else:
        server_key = server[server.find('/')+2:server.find('.')]# beta: "http://v25rc2.demo.encodedcc.org"
    if server_key != 'www':
        acc_key = server_key + '_accession'
    return acc_key
    
def file_set_property(fid,key,value,proj_id=None,add_only=False,test=False,verbose=False):
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
                print >> sys.stderr, "Note: file %s already has property '%s' set to '%s'" % (path,key,properties[key])
            return properties[key]
        elif add_only:
            if verbose:
                print >> sys.stderr, "Error: file %s already has property '%s' and is not being updated from '%s'" % (path,key,properties[key])
            return properties[key]
        elif verbose:
                print >> sys.stderr, "Warning: file %s has property '%s' but is being updated to '%s'" % (path,key,value)
    properties[key] = value
    if test:
        if verbose:
            print >> sys.stderr, "  - Test set %s with %s='%s'" % (path,key,value)
    else:
        dxfile.set_properties(properties)
        if verbose:
            print >> sys.stderr, "  - set %s with %s='%s'" % (path,key,value)
    return properties[key]
    
def umbrella_folder(folder,default,proj_name=None,exp_type=None,sub_folder=None,genome=None,annotation=None):
    '''Returns a normalized umbrella folder (that holds the experiments of a given type).'''
    if folder != default:
        return folder_normalize(folder)
    if sub_folder != None:
        sub_folder = folder_normalize(sub_folder,starting=False)
            
    # No change to default, so build from parts if available
    if exp_type == None:
        return folder_normalize(folder)

    if proj_name == PRODUCTION_PROJECT:
        if sub_folder == None:
            sub_folder = "runs/" 
        if exp_type == "long-rna-seq":
            folder = "/long-RNA-seq/"
        elif exp_type == "small-rna-seq":
            folder = "/small-RNA-seq/"
        elif exp_type == "dnase-seq":
            folder = "/DNAse-seq/runs/"
        elif exp_type == "dna-me":
            folder = "/WG Bisulfite (Methylation)/"
        elif exp_type == "chip-seq":
            folder = "/ChIP-seq/"
        else:
            folder = "/" + exp_type + '/'
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

    if sub_folder != None:
        folder += sub_folder

    if genome != None:
        folder +=  genome + '/'
        if annotation != None and genome == 'mm10':
            folder += annotation + '/'
            
    return folder


def duration_string(total_seconds,include_seconds=True):
    '''Returns a duration string.'''
    try: 
        secs = int(total_seconds)
    except:
        return None
    m, s = divmod(secs, 60)
    h, m = divmod(m, 60)
    d, h = divmod(h, 24)
    y, d = divmod(d, 365)
    
    if include_seconds:
        if y > 0:
            return "%dy%03dd%02dh%02dm%02s" % (y, d, h, m, s)
        elif d > 0:
            return      "%dd%02dh%02dm%02s" % (d, h, m, s)
        elif h > 0: 
            return          "%dh%02dm%02ds" % (h, m, s) 
        elif m > 0: 
            return               "%dm%02ds" % (m, s)
        else: 
            return                    "%ds" % (s)
    else:
        if s >= 30:
            m +=1
        if y > 0:
            return "%dy%03dd%02dh%02dm" % (y, d, h, m)
        elif d > 0:
            return      "%dd%02dh%02dm" % (d, h, m)
        elif h > 0: 
            return           "%dh%02dm" % (h, m) 
        else: 
            return                "%dm" % (m)


def format_duration(beg_seconds,end_seconds,include_seconds=True):
    '''Returns formatted string difference between two times in seconds.'''
    duration = end_seconds - beg_seconds
    return duration_string(duration,include_seconds)

