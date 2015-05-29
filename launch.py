#!/usr/bin/env python2.7
# launch.py 0.0.1

import argparse,os, sys, json
#import urlparse, subprocess, itertools, logging
from datetime import datetime
from collections import deque

import dxpy
import dxencode

### TODO:
#   1) Launchers should create the workflow_run encoded objects.  REQUIRED new status: 'launch_error', 'started', 'failed', 'partial_success')
#   2) Launchers should be able to support 3 level combinations [br:[tr,tr],br:[tr,tr]] -> [br,br] -> cr (DNase needs)
#   3) Launchers may fail, so how to report this to encoded?

# NOTES: This command-line utility will run the a pipeline for a single replicate or combined reps.
#      - All results will be written to a folder /<results_folder>/<exp_id>/rep<#>_<#>.
#      - If any results are already in that directory, then the steps that created those results
#        will not be rerun.
#      - If any jobs for the experiment and replicate are already running, nothing new will be
#        launched.
#      - Almost all the code is generic and derived classes will only need to entend or override a 
#        a few things.  It relies upon JSON pipeline descriptions which in simple cases can be 
#        directly read from dx (e.g. srnaLaunch.py).  More compex examples (e.g. lrnaLaunch.py)
#        must hard-code the pipeline json to handle repeated apps and dx name conflicts.
#        Tokens are used to abstract dx.app input/outout file names to avoid collisions.
#        - STEP_ORDER is the list of steps in the pipeline
#        - REP_STEPS and COMBINED_STEPS contain step definitions and enforces dependencies by input
#          file tokens matching to result file tokens of earlier steps.
#        - FILE_GLOBS is needed for locating result files from prior runs.
#      - Combined replicate processing is optional and depends upon COMBINED_STEP, etc.
#        (e.g. rampageLaunch.py).
#      - Control files are supported but most likely will require function overrides in derived
#        classes (e.g. rampageLaunch.py)
#      - In REP_STEPS and COMBINED_STEPS, both 'inputs' and 'results' are key:value lists where the key is used
#        to attach results of one step to the inputs of another, while values must match dxapp.json names.
#        Further, any input key ending in '_set' is expected to be an array of files (e.g. "input": { "reads_set": "reads" }  
#        will expect one or more files to be found and used for dxapp.json input named "reads" ).
#      - The FILE_GLOBS are matched by input/result keys, not values.
#        (e.g. "input": { "pooled_bam": "bam" } matches file_glob "pooled_bam" and dxapp.json input named "bam".
#        This example would allow the same applet to be used on replicate bams and poold bams in the same pipeline.)
#
# Example derived classes: long-rna-seq-pipeline: lrnaLaunch.py, srnaLaunch.py rampageLaunch.py

class Launch(object):
    '''
    Launch module is the common code for pipeline specific launcher classes.  
    Launchers dynamically create and launch a workflow covering any needed steps of a pipeline.
    Steps needed are thos that have available inputs but not already generatedresults.
    '''

    PIPELINE_NAME = "MUST REPLACE IN DREIVED CLASS"
    ''' This must match the assay type returned by dxencode.get_assay_type({exp_id}).'''

    PIPELINE_HELP = "Launches '"+PIPELINE_NAME+"' pipeline " + \
                    "analysis for one replicate or combined replicates. "
    ''' This help title should name pipline and whether combined replicates are supported.'''
                    
    SERVER_DEFAULT = 'www'
    '''At this time there is no need to use the any but the one true server for launching.'''
    
    FOLDER_DEFAULT = '/runs/'
    ''' This the default location to place results folders for each experiment.'''
    
    REP_STEP_ORDER = None # [ "MUST","REPLACE","IN","DREIVED","CLASS" ]
    '''The (artifically) linear order of replicate-only pipeline steps.'''
    COMBINED_STEP_ORDER = None  # [ "MUST","REPLACE", "IF", "COMBINED", "STEPS", "EXIST" ]
    '''The (artifically) linear order of combined replicate pipeline steps.'''

    REP_STEPS = None # { "MAY": {}, "NEED": {}, "TO": {}, "ADD": {}, "IF": "NOT SIMPLE" }
    '''Replicate step objects which can be discoverd in simple piplines.'''
    COMBINED_STEPS = None # { "NEED": {}, "TO": {}, "REPLACE": {}, "IF": "COMBINED IS SUPPORTED" }
    '''Combined step objects which are needed if the pipeline has combined replicate processing.'''
    # STEP objects contain the following:
    #    "step-name":   { # Must match STEP_ODER and may match dx applet
    #        "app":     "dx-applet-name", # Must match applet
    #        "params":  { "param_tag": "dx_input_name", ...}, # Non-files: Tag must match psv key.
    #        "inputs":  { "input_tag": "dx_input_name", ...},  # Files: Tag *may match* name
    #        "results": { "result_tag": "dx_output_name", ...}} # Files: Tag *may match* name
    # Note, tags *may* be the same as dx names, but MUST MATCH other tags to define dependencies as: 
    #     stepA:result_tag == stepB:input_tag
    # Also result_tags must be in FILE_GLOB.keys()
    #
    # NOTE: Simple pipelines have no reused steps or dx name to param complications and the 
    #       STEP objects can be generated directly from dx apps.  If unsure, print the json from
    #       self.build_simple_steps(verbose=True) to determine if it will work for your pipeline.
    
    FILE_GLOBS = {"MAY": {}, "NEED": {}, "TO": {}, "REPLACE": {} }
    '''Globs for discovering existing results. Can be discoverd in simple piplines.'''

    CONTROL_FILE_GLOB = None
    '''Glob for discovering control files if experiment has controls. Othewise: None'''

    REFERENCE_FILES = {
        # For looking up reference file names.
        }

    GENOMES_SUPPORTED = ['hg19', 'mm10']
    GENOME_DEFAULT = 'hg19'
    ''' This the default Genome that long RNA-seq experiments are mapped to.'''
    
    STANDARD_LABELS = {
        'long-rna-seq':  { 'long': 'Long-RNA-seq',  'short': 'lrna',  'se_and_pe': True  }, 
        'small-rna-seq': { 'long': 'Small-RNA-seq', 'short': 'srna',  'se_and_pe': False }, 
        'rampage':       { 'long': 'Rampage',       'short': 'ramp',  'se_and_pe': False }, 
        'dnase-seq':     { 'long': 'DNase-seq',     'short': 'dnase', 'se_and_pe': True  },
    }
    '''Standard labelling requires exp_type specific labels.  This can be overridden in descendent classes.'''

    def __init__(self):
        '''
        Launch expects to be the base class for derived pipeline specific launch classes.
        As such it will not run stand alone.
        an experiment.
        '''
        self.args = {} # run time arguments
        self.server_key = self.SERVER_DEFAULT
        self.proj_name = None
        self.project = None
        self.proj_id = None
        self.exp = {}  # Will hold the encoded exp json
        self.psv = {} # will hold pipeline specific variables.
        self.multi_rep = False       # This run includes more than one replicate
        self.combined_reps = False   # This run includes combined replicate workflows
        print # TEMPORARY: adds a newline to "while retrieving session configuration" unknown error
    
    def get_args(self,parse=True):
        '''Parse the input arguments.'''
        # MAY NEED TO EXTEND OR REPLACE in derived class
        # see long-rna-seq: lrnaLaunch.py for and example of extending
        
        # Change this description to PIPELINE SPECIFIC version
        ap = argparse.ArgumentParser(description=self.PIPELINE_HELP +
                    "Can be run repeatedly and will launch only the steps that are needed to " +
                    "finish the pipeline. All results will be placed in the folder: " +
                    "/<resultsLoc>/<experiment>/<replicate>.")

        ap.add_argument('-e', '--experiment',
                        help='ENCODED experiment accession',
                        required=True)

        ap.add_argument('-br','--br', '--biological-replicate',
                        help="Biological replicate number",
                        type=int,
                        #default='1',
                        required=False)

        ap.add_argument('-tr','--tr', '--technical-replicate',
                        help="Technical replicate number (default: 1)",
                        type=int,
                        default='1',
                        required=False)

        if self.CONTROL_FILE_GLOB != None:
            # Include this argument for pipelines with controls
            ap.add_argument('-c', '--control',
                            help='The control file for this experiment.',
                            required=False)

        if self.COMBINED_STEPS != None:
            # Include this argument for pipelines with combined replicate steps
            ap.add_argument('-cr','--cr','--combine-replicates',
                            help="Combine or compare two replicates (e.g.'1 2_2').'",
                            nargs='+',
                            required=False)

        ap.add_argument('--project',
                        help="Project to run analysis in (default: '" + \
                                                        dxencode.env_get_current_project() + "')",
                        required=False)

        ap.add_argument('--refLoc',
                        help="The location to find reference files (default: '" + \
                            dxencode.REF_PROJECT_DEFAULT + ":" + dxencode.REF_FOLDER_DEFAULT + "')",
                        default=dxencode.REF_FOLDER_DEFAULT,
                        required=False)

        ap.add_argument('-f','--folder',
                        help="The location to to place results folders (default: '<project>:" + \
                                                                  self.FOLDER_DEFAULT + "')",
                        default=self.FOLDER_DEFAULT,
                        required=False)

        ap.add_argument('--run',
                        help='Run the workflow after assembling it.',
                        action='store_true',
                        required=False)

        ap.add_argument('--server',
                        help="Server to post files to (default: '" + self.SERVER_DEFAULT + "')",
                        default=self.SERVER_DEFAULT,
                        required=False)

        ap.add_argument('--test',
                        help='Test run only, do not launch anything.',
                        action='store_true',
                        required=False)

        ap.add_argument('--force',
                        help='Force rerunning all steps.',
                        action='store_true',
                        required=False)

        # Not yet supported.
        #ap.add_argument('-x', '--export',
        #                help='Export generic Workflow (no inputs) to DNA Nexus project',
        #                action='store_true',
        #                required=False)

        if parse:
            return ap.parse_args()
        else:
            return ap

    def pipeline_specific_vars(self,args,verbose=False):
        '''Adds pipeline specific variables to a dict, for use building the workflow.'''
        # PLEASE EXTEND in derived class

        # psv can contain any variables, but it must contain these at a minimum:
        # - Any non-file input param needed to launch the workflow
        # - 'resultFolder' - full dx path (without project) to the results folder of the specific run
        # - 'name' - A short name used for specific workflow run. 
        # - 'description' - generic description of the pipeline.
        # - 'title'/['subtitle'] for command line output announcing what will be done 

        # Start with dict containing common variables
        print "Retrieving experiment specifics..."
        self.server_key = args.server
        psv = self.common_variables(args,key=self.server_key)
        if psv['exp_type'] != self.PIPELINE_NAME:
            print "Experiment %s is not for '%s' but for '%s'" \
                                           % (psv['experiment'],self.PIPELINE_NAME,psv['exp_type'])
            sys.exit(1)
        
        # Now add pipline specific variables here...


        if verbose:
            print "Pipeline Specific Vars:"
            print json.dumps(psv,indent=4)
        return psv


    def find_ref_files(self,priors):
        '''Locates all reference files based upon organism and gender.'''
        # PLEASE REPLACE in derived class
        assert "must replace!" and False
        
        # Add your pipeline specific referernce file detection here
        
        # Examples:
        # If all the ref files are in json keyed by ['genome']['gender'] then try:
        for ref in self.REFERENCES_FILES.keys():
            dxfile = self.psv['refLoc']+self.REFERENCES_FILES[ref][self.psv['genome']][self.psv['gender']]
            fid = dxencode.find_file(dxfile,dxencode.REF_PROJECT_DEFAULT)
            if fid == None:
                sys.exit("ERROR: Unable to locate ref file: '" + dxfile + "'")
            else:
                priors[ref] = fid
        
        # Or explicitly by key
        chrom_sizes = self.psv['refLoc']+self.REFERENCE_FILES['chrom_sizes'][self.psv['genome']][self.psv['gender']]
        chrom_sizes_fid = dxencode.find_file(chrom_sizes,dxencode.REF_PROJECT_DEFAULT)
        if chrom_sizes_fid == None:
            sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chrom_sizes + "'")
        else:
            priors['chrom_sizes'] = chrom_sizes_fid
        self.psv['ref_files'] = self.REFERENCE_FILES.keys()
        

    def find_control_file(self,rep,default=None):
        '''Attempts to find an appropriate control file.'''
        # NEED TO REPLACE for pipelines with controls
        # Find example in long-rna-seq-pipeline: rampageLaunch.py

        if self.CONTROL_FILE_GLOB == None:
            return None
            
        
    def find_all_control_files(self,test):
        '''Attempts to find all control files needed for a run
        .'''
        # MAY NEED TO REPLACE for pipelines with controls

        if self.CONTROL_FILE_GLOB == None:
            return
        print "Checking for control files..."
        default_control = None
        if 'control' in self.psv:
            default_control = self.psv['control'] # akward but only rep1 may be run
        for rep in self.psv['reps'].values():
            control = self.find_control_file(rep,default_control)
            rep['inputs']['Control'] = dxencode.find_and_copy_read_files(rep['priors'], \
                                                    [ control ], test, 'control_bam', \
                                                    rep['resultsFolder'], False, self.proj_id)

    def find_combined_inputs(self,steps,file_globs):
        '''Finds the inputs for a combined run in the directories of the replicate runs.'''
        # MAY NEED TO REPLACE for pipelines with combined steps
        
        if not self.combined_reps:
            return

        print "Checking for combined-replicate inputs..."
        inputs = {}
        for step in self.psv['path']:
            for file_token in steps[step]['inputs'].keys():
                if file_token[-2] != '_':
                    continue
                rep_key = file_token[-1].lower()
                if rep_key not in self.psv['reps'].keys():
                    print "*** " + rep_key + " not found in psv['reps']"
                    continue
                rep = self.psv['reps'][rep_key]
                # TODO: No need for multiples at this time.  Deal with it when it comes up.
                fid = dxencode.find_file(rep['resultsFolder'] + file_globs[file_token],\
                                                    self.proj_id, multiple=False, recurse=False)
                if fid != None:
                    self.psv['priors'][file_token] = fid
                    inputs[file_token] = [ fid ]
                elif not self.multi_rep:
                    print "Error: Necessary '%s' for combined run, not found in '%s'." \
                                                               % (file_token, rep['resultsFolder'])
                    print "       Please run for single replicate first."
                    sys.exit(1)
        
        self.psv['inputs'] = inputs
    

    ############## NOT EXPECTED TO OVERIDE THE FOLLOWING METHODS ############
    def common_variables(self,args,fastqs=True,key='default'):
        '''Initializes dict with common variables from args and encoded.'''
        # NOT EXPECTED TO OVERRIDE
        
        cv = {}
        self.proj_name = dxencode.env_get_current_project()
        if self.proj_name == None or args.project != None:
            self.proj_name = args.project
        if self.proj_name == None:
            print "Please enter a '--project' to run in."
            sys.exit(1)
        self.project = dxencode.get_project(self.proj_name)
        self.proj_id = self.project.get_id()        

        cv['project']    = self.proj_name
        cv['experiment'] = args.experiment
        
        ### TODO: get all replicates, not relying on -cr, -br, -tr
        self.exp = dxencode.get_exp(cv['experiment'],key=key)
        cv['exp_type'] = dxencode.get_assay_type(cv['experiment'],self.exp)
        full_mapping = dxencode.get_full_mapping(cv['experiment'],self.exp)
        cv['reps'] = self.load_reps(args, cv, cv['experiment'], self.exp, full_mapping)
        
        assert 'a' in cv['reps']

        # Only supported genomes
        cv['gender'] = cv['reps']['a']['sex']
        organism = cv['reps']['a']['organism']
        if organism in dxencode.GENOME_DEFAULTS:
            cv['genome'] = dxencode.GENOME_DEFAULTS[organism]
        else:
            print "Organism %s not currently supported" % organism
            sys.exit(1)

        # Paired ends?
        if self.combined_reps and cv['reps']['a']['paired_end'] != cv['reps']['b']['paired_end']:
            print "Replicates are expected to be both paired-end or single-end!  Check encoded."
            sys.exit(1)
        cv['paired_end'] = cv['reps']['a']['paired_end']

        # Default locations
        cv['refLoc'] = args.refLoc
        if cv['refLoc'] == dxencode.REF_FOLDER_DEFAULT:
            cv['refLoc'] = dxencode.REF_FOLDER_DEFAULT + cv['genome'] + '/'
        if not cv['refLoc'].endswith('/'):
            cv['refLoc'] += '/' 
        cv['resultsLoc'] = dxencode.umbrella_folder(args.folder,self.FOLDER_DEFAULT,self.proj_name,cv['exp_type'],cv['genome'])
        cv['resultsFolder'] = cv['resultsLoc'] + cv['experiment'] + '/'
        for ltr in cv['reps'].keys():
            cv['reps'][ltr]['resultsFolder'] = cv['resultsLoc'] + cv['experiment'] + '/' + cv['reps'][ltr]['rep_tech'] + '/'
        
        # Standard labels:
        if cv['exp_type'] in self.STANDARD_LABELS:
            labels = self.STANDARD_LABELS[cv['exp_type']]
            cv['description'] = "The ENCODE "+labels['long']+" pipeline"
            cv['title'] = labels['long']
            cv['name'] = labels['short'] + "_"+cv['genome']
            if cv['gender'] == 'female':
                cv['name'] += "XX"
            else:
                cv['name'] += "XY"
            if labels['se_and_pe']:
                if cv['paired_end']:
                    cv['title'] += " paired-end"
                    cv['name'] += "PE"
                else:
                    cv['title'] += " single-end"
                    cv['name']  += "SE"
            cv['title'] += " " + cv['experiment']
            cv['name']  += "_" + cv['experiment']
            for ltr in cv['reps'].keys():
                rep = cv['reps'][ltr]
                rep['title'] = cv['title'] + " - "+rep['rep_tech']
                if 'library_id' in rep:
                    rep['title'] += " (library '"+rep['library_id']+"')"
                rep['title'] += " on "+cv['genome']+", "+cv['gender']
                rep['name']  = cv['name']  + rep['rep_tech']
            cv['title']   += " - "+cv['rep_tech']+" on "+cv['genome']+", "+cv['gender']
            cv['name']    += "_"+cv['rep_tech']

        return cv

    def load_reps(self,args, cv, exp_id, exp=None, full_mapping=None):
        '''
        Gathers replicates from encoded, then builds the cv['reps'] list from requested or all.
        '''
        # NOT EXPECTED TO OVERRIDE
        cv_reps = {}
        reps = dxencode.get_reps_from_enc(exp_id, load_reads=True, exp=exp, full_mapping=full_mapping)
        rep_techs = []
        if 'cr' in args and args.cr != None:
            self.multi_rep = True
            if len(args.cr) != 2:
                print "Specify which two replicates to compare (e.g. '1_2 2')."
                sys.exit(1)
            self.combined_reps = True
            for rep_tech in args.cr:
                # Normalize rep_tech string
                if '_' not in rep_tech:
                    rep_tech = rep_tech + '_1' # default to tech_rep 1
                if not rep_tech.startswith('rep'):
                    rep_tech = 'rep' + rep_tech
                rep_techs.append(rep_tech)
        elif args.br != None and args.br != 0:
            rep_tech = 'rep' + str(args.br) + '_' + str(args.tr)
            rep_techs.append(rep_tech)
            # TODO: support arrays for br, tr
            #self.combined_reps = False
            #for ix,br in enumerate(args.br):
            #    tr = 1
            #    if ix < len(args.tr):
            #        tr = args.tr[ix]
            #    rep_tech = 'rep' + str(br) + '_' + str(tr)
            #    rep_techs.append(rep_tech)
            
        # Find each requested rep_tech in reps
        letters = deque('abcdefghijklmnopqrstuvwxyz')
        if len(rep_techs) > 0:
            # Add specifically requested rep_techs to lettered slots
            self.multi_rep = (len(rep_techs) > 1)
            for rep_tech in sorted( rep_techs ):
                found_ix = -1
                for ix,rep in enumerate(reps):
                    if rep['rep_tech'] == rep_tech:
                        found_ix = ix
                if found_ix == -1:
                    for ltr in cv_reps.keys():
                        if cv_reps[ltr]['rep_tech'] == rep_tech:
                            print "Specify different replicates to compare (e.g. 'rep1_1 rep2_1')."
                            sys.exit(1)
                    print "ERROR: requested '"+rep_tech+"' not found in ENCODE " + cv['experiment']
                    sys.exit(1)
                ltr = letters.popleft()
                cv_reps[ltr] = reps[found_ix]
                del reps[found_ix] # rep is only used once!
        else:
            # Assume combined if exactly 2 replicates are available
            for rep in reps:
                ltr = letters.popleft()
                cv_reps[ltr] = rep
            # Assume combined reps if supported AND exactly 2 reps AND for different biological reps
            if self.COMBINED_STEP_ORDER != None and len( self.COMBINED_STEP_ORDER ) > 0:
                if len(reps) == 2 and cv_reps['a']['br'] != cv_reps['b']['br']:
                    self.combined_reps = True
            self.multi_rep = (len(reps) > 1)
            
        # NOTE: Plans for combining more than 2 reps and multi-levels
        # Currently reps are keyed on single letter (e.g. 'a') and added to workflow in sort order
        # Combined level keying should build workflows in order AND allow tying results to inputs of combined steps
        # Combined reps need to be added after last rep in combination, so just add to the end:
        # Example reps['a'] and reps['b'] are combined in reps['z-ab']
        # However, more complex levels still could need to be combined it is best to have combined reps point
        # their dependent reps.  reps['z-ab']['parents'] = ['a','b']
        # Notice that psv is essentially reps['z-ab']
            
        # mult-rep rep_tech: 
        if self.combined_reps:
            cv['rep_tech'] = 'reps' + cv_reps['a']['rep_tech'][3:] + \
                                '-' + cv_reps['b']['rep_tech'][3:]  # combined delimited by '-'
            cv['tributaries'] = ['a','b']
        elif self.multi_rep:
            cv['rep_tech'] = 'reps:' + cv_reps['a']['rep_tech'][3:]
            for ltr in sorted(cv_reps.keys()):
                if ltr != 'a':
                    cv['rep_tech'] += ','+cv_reps[ltr]['rep_tech'][3:] # not combined are delimited by ','
        else:
            cv['rep_tech'] = cv_reps['a']['rep_tech'] 

        # A little more rep tidying
        for ltr in cv_reps.keys():
            rep = cv_reps[ltr]
            rep['concat_id'] = 'reads'
            if rep['paired_end']: 
                rep['concat_id2'] = 'reads2'
                
        return cv_reps


    def build_a_step(self, applet, file_globs, proj_id):
        ''' create input object for a step and extends the file_globs dict as appropriate.'''
        # NOT EXPECTED TO OVERRIDE
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


    def build_simple_steps(self,pipe_path, proj_id, verbose=False):
        '''
        Builds dict of steps for the apps in the pipeline and a dict of file_globs for look up.
        Only works for pipelines where every step in pipe_path is a distinct app,
        and each result glob is uniq for all pipeline results.
        '''
        # NOT EXPECTED TO OVERRIDE
        steps = {}
        file_globs = {}
        for step in pipe_path:
            steps[step] = self.build_a_step(step, file_globs, proj_id)
        if verbose:
            print "STEPS = "
            print json.dumps(steps,indent=4)
            print "FILE_GLOBS = "
            print json.dumps(file_globs,indent=4)

        return [ steps, file_globs ]


    def assemble_steps_and_globs(self,pipe_path,combined_steps=False):
        '''Assemble the requested steps and file globs.'''
        # NOT EXPECTED TO OVERRIDE
            
        if self.FILE_GLOBS != None and self.REP_STEPS != None:
            # Use for complex pipelines:
            file_globs = self.FILE_GLOBS
            if not combined_steps:
                pipe_steps = self.REP_STEPS
                if type(pipe_path) == dict:
                    assert 'se' in pipe_path and 'pe' in pipe_path
                    if self.psv['paired_end']:
                        pipe_path = pipe_path['pe']
                    else:
                        pipe_path = pipe_path['se']
            elif self.COMBINED_STEPS:
                pipe_steps = self.COMBINED_STEPS
            else:
                print "ERROR: Could not find COMBINED_STEPS, though combined replicated requested."
                sys.exit(1)
        else:
            # Use for simple pipelines
            pipe_steps, file_globs = self.build_simple_steps(pipe_path,self.proj_id)

        # psv must reference pipe_path
        if combined_steps:
            self.psv['path'] = pipe_path
        else:
            for rep in self.psv['reps'].values():
                rep['path'] = pipe_path
            
        return [ pipe_steps, file_globs ]

    def find_prior_results(self,pipe_path,steps,results_folder,file_globs):
        '''Looks for all result files in the results folder.'''
        priors = {}
        for step in pipe_path:
            for file_token in steps[step]['results'].keys():
                fid = dxencode.find_file(results_folder + file_globs[file_token],self.proj_id, \
                                                                                    recurse=False)
                if fid != None:
                    priors[file_token] = fid
        return priors


    def finding_rep_inputs_and_priors(self,steps,globs,test):
        '''Finds the inputs and priors for a run.'''
        # NOT EXPECTED TO OVERRIDE
        
        print "Checking for prior results..."
        # NOTE: priors is a dictionary of fileIds that will be used to determine stepsToDo
        #       and fill in inputs to workflow steps
        for rep in self.psv['reps'].values():
            if not test:
                if not dxencode.project_has_folder(self.project, rep['resultsFolder']):
                    self.project.new_folder(rep['resultsFolder'],parents=True)
            rep['priors'] = self.find_prior_results(rep['path'],steps,rep['resultsFolder'],globs)

        print "Checking for input files..."
        # Find all reads files and move into place
        # TODO: files could be in: dx (usual), remote (url e.g.https://www.encodeproject.org/...)
        #       or possibly local, Currently only DX locations are supported.
        for rep in self.psv['reps'].values():
            rep['inputs'] = {}
            reads_token = 'reads'
            if self.psv['paired_end']:
                reads_token += '1'
            rep['inputs']['Reads1'] = dxencode.find_and_copy_read_files(rep['priors'], \
                                                rep['fastqs']['1'],  test, reads_token, \
                                                rep['resultsFolder'], False, self.proj_id)
            # Note: rep['fastqs']['2'] and rep['inputs']['Reads2'] will be empty on single-end
            rep['inputs']['Reads2'] = dxencode.find_and_copy_read_files(rep['priors'], \
                                                rep['fastqs']['2'], test, 'reads2', \
                                                rep['resultsFolder'], False, self.proj_id)


    def finding_combined_inputs_and_priors(self,steps,globs,test):
        '''Finds the inputs and priors for a combined run.'''
        # NOT EXPECTED TO OVERRIDE
        
        if not self.combined_reps:
            return
        
        print "Checking for combined-replicate priors..."
        if not test:
            # should be a noop as replicate level folders were already created
            if not dxencode.project_has_folder(self.project, self.psv['resultsFolder']):
                self.project.new_folder(self.psv['resultsFolder'],parents=True)
        self.psv['priors'] = self.find_prior_results(self.psv['path'],steps,self.psv['resultsFolder'],globs)

        # Checking for combined input files...
        # The following step may be overwritten in descendent class
        self.find_combined_inputs(steps,globs)
          
            
    def find_all_ref_files(self):
        '''Locates all reference files based upon organism and gender.'''
        # NOT EXPECTED TO OVERRIDE
        
        print "Looking for reference files..."
        for rep in self.psv['reps'].values():
            # Need multiple copies of control files because they are put into priors!
            self.find_ref_files(rep['priors']) # pipeline specific
        if self.combined_reps:
            self.find_ref_files(self.psv['priors'])


    def determine_steps_to_run(self,pipe_path, steps, priors, deprecate, force=False, verbose=False):
        '''Determine what steps need to be done, base upon prior results.'''
        # NOT EXPECTED TO OVERRIDE
        will_create = []
        steps_to_run = []
        for step in pipe_path:
            # Force will include the first step with all its inputs
            # This should avoid forcing concat if it isn't needed
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
                                                        project=self.proj_id, return_handler=False)
            if dxApp == None:
                print "ERROR: failure to locate app '"+app+"'!"
                sys.exit(1)

        return steps_to_run


    def determine_steps_needed(self,replicate_steps, combined_steps, force=False):
        '''Determine steps needed for replicate(s) and combined, base upon prior results.'''
        # NOT EXPECTED TO OVERRIDE
        
        print "Determining steps to run..."
        # NOTE: stepsToDo is an ordered list of steps that need to be run
        for rep in self.psv['reps'].values():
            rep['steps'] = replicate_steps
            rep['deprecate'] = [] # old results will need to be moved/removed if step is rerun
            rep['stepsToDo'] = self.determine_steps_to_run(rep['path'], replicate_steps, \
                                                    rep['priors'], rep['deprecate'], force=force)
        if self.combined_reps:
            self.psv['steps'] = combined_steps
            self.psv['deprecate'] = [] # old results will need to be moved/removed if step is rerun
            self.psv['stepsToDo'] = self.determine_steps_to_run(self.psv['path'], combined_steps, \
                                            self.psv['priors'], self.psv['deprecate'], force=force)


    def find_resuts_from_prev_branch(self, this_branch, branches, inp_token, expect_set):
        '''
        If combined_reps step, look for input from results of replicate level step.
        '''
        if not self.combined_reps:    # Not even doing a wf with combining branches 
            return None
        if 'tributaries' not in this_branch:  # Not a combining branch
            return None
        if inp_token not in self.FILE_GLOBS: # assertable
            return None
        # The combined_step input tokens should either end in '_a', '_b' or the like,
        # Or the step expects to combine N files so expect_set=True and all matching results should be combined.
        
        # Using letters since 3 tech_rep branches could be combined into 2 bio_rep branches and then 1 experiment level branch
        # Thus tr branch_ids would be: 'a','b','c','d','e','f'
        #  and br branch_ids would be 'za','zb', and those branches flow into psv (or 'zc').
        # So if br branch 'za' expects 3 inputs ending '_a','_b','_c' then
        #       br branch 'zb' would expect the same 3, not '_d','_e','_f'
        letters = deque('abcdefghijklmnopqrstuvwxyz')
        looking_for = None
        if expect_set:
            # TODO: Establish expectation that joingin results either end in '_a', etc. or '_set'
            # TODO: Refactor naming to refer to reps, branches and rarely "runs' 
            if not inp_token.endswith('_set'):
                return None
            #looking_for = 'set'  # How to handle sets?
            results_array = []
        elif inp_token[-2] == '_' and inp_token.lower()[-1] in letters:
            looking_for = inp_token.lower()[-1]
        if looking_for == '': # Can we do this
            return None
            
        for branch_id in this_branch['tributaries']:
            ltr = letters.popleft()
            if looking_for != None and looking_for != ltr:
                continue
            if branch_id not in branches:  # assertable
                continue
            if 'prevStepResults' not in branches[branch_id]:
                continue
            prev_results = branches[branch_id]['prevStepResults']
            # How to match prevResults[file_token='bwa_bam'] in prevResults to compare to step['inputs'][file_token='bam_a']
            # The only way is by matching file globs which MUST tie branch leaps together.
            for result_token in prev_results.keys():
                if result_token in self.FILE_GLOBS:
                    if self.FILE_GLOBS[result_token] == self.FILE_GLOBS[inp_token]:
                        if expect_set:
                            results_array.append(prev_results[result_token])
                            break  # Done with this branch, but look for more results in next branches
                        if looking_for == ltr:
                            return prev_results[result_token] # found result on branch matched by letter
                        else:
                            return prev_results[result_token] # return first match from first branch
        
        # Went through all branches so return any results found
        if not expect_set or len(results_array) == 0:
            return None
        return results_array


    def create_or_extend_workflow(self,run, run_id, wf=None,app_proj_id=None,test=False,template=False):
        '''
        This function will populate a workflow for the steps in run['stepsToDo'] and return 
        the workflow unlaunched.   It relies on steps dict which contains input and output
        requirements, pvs (pipeline specific variables) dictionary and run (run specific
        dictionary), which contains input and previous results already in results dir
        '''
        # NOT EXPECTED TO OVERRIDE

        if len(run['stepsToDo']) < 1:
            return None
        if app_proj_id == None:
            app_proj_id = self.proj_id

        # create a workflow object
        if self.multi_rep:
            wf_folder = self.psv['resultsFolder']
            wf_name = self.psv['name']
        else:
            wf_name = run['name']
            wf_folder = run['resultsFolder']
        if not test and wf == None:
            wf = dxpy.new_dxworkflow(title=wf_name,name=wf_name,folder=wf_folder,
                                           project=self.proj_id,description=self.psv['description'])

        # NOTE: prevStepResults dict contains links to result files to be generated by previous steps
        if 'prevStepResults' not in run:
            run['prevStepResults'] = {}
        steps = run['steps']
        for step in run['stepsToDo']:
            appName = steps[step]['app']
            app = dxencode.find_applet_by_name(appName, app_proj_id)
            inp_defs = app.describe().get('inputSpec') or []
            appInputs = {}
            # file inputs
            for fileToken in steps[step]['inputs'].keys():
                appInp = steps[step]['inputs'][fileToken]
                if fileToken in run['prevStepResults']:
                    appInputs[ appInp ] = run['prevStepResults'][fileToken]
                else:
                    for inp_def in inp_defs:
                        if inp_def["name"] == appInp:
                            break
                    assert(inp_def["name"] == appInp)
                    # Check to see if results from a previous branch (e.g. rep level pipeline) is input to the combining step
                    expect_set = inp_def["class"] == 'array:file'
                    prev_results = self.find_resuts_from_prev_branch(run, self.psv['reps'], fileToken, expect_set)
                    if prev_results != None:
                        appInputs[ appInp ] = prev_results
                    elif fileToken in run['priors']:
                        if inp_def["class"] == 'array:file':
                            if isinstance(run['priors'][fileToken], list):
                                appInputs[ appInp ] = []
                                for fid in run['priors'][fileToken]:
                                    appInputs[ appInp ] += [ dxencode.FILES[fid] ]
                            else:
                                appInputs[ appInp ] = [ dxencode.FILES[ run['priors'][fileToken] ] ]
                        else:
                            assert(not isinstance(run['priors'][fileToken], list))
                            appInputs[ appInp ] = dxencode.FILES[ run['priors'][fileToken] ]
                if appInp not in appInputs and not template:
                    print "ERROR: step '"+step+"' can't find input '"+fileToken+"'!"
                    sys.exit(1)
            # Non-file app inputs
            if 'params' in steps[step]:
                for param in steps[step]['params'].keys():
                    appParam = steps[step]['params'][param]
                    if param in run:
                        appInputs[ appParam ] = run[param]
                    elif param in self.psv:
                        appInputs[ appParam ] = self.psv[param]
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
                    run['prevStepResults'][ fileToken ] = 'fake-for-testing'
                else:
                    run['prevStepResults'][ fileToken ] = dxpy.dxlink({'stage': stageId, \
                                                                'outputField': appOut })

        if test:
            return None
        else:
            return wf


    def report_run_plans(self,run=None):
        '''Report the plans before executing them.'''
        # NOT EXPECTED TO OVERRIDE

        if run == None:
            run = self.psv
        print "Running: "+run['title']
        if 'subTitle' in run:
            print "         "+run['subTitle']

        # Inputs:
        if self.multi_rep:
            for ltr in sorted( self.psv['reps'].keys() ):
                rep = self.psv['reps'][ltr]
                for input_type in sorted( rep['inputs'].keys() ):
                    if len(rep['inputs'][input_type]) > 0: 
                        print "- " + input_type + " " + rep['rep_tech'] +':'
                        for fid in rep['inputs'][input_type]:
                            print "  " + dxencode.file_path_from_fid(fid)
        # single-rep and combined work both do run:
        if not self.multi_rep or self.combined_reps:
            for input_type in sorted( run['inputs'].keys() ):
                if len(run['inputs'][input_type]) > 0: 
                    print "- " + input_type + " " + run['rep_tech'] +':'
                    for fid in run['inputs'][input_type]:
                        print "  " + dxencode.file_path_from_fid(fid)

        print "- Reference files:"
        # Should be enough to show the ref priors from 'a' rep for single-rep, multi or combined
        for token in self.psv['ref_files']:
            print "  " + dxencode.file_path_from_fid(self.psv['reps']['a']['priors'][token],True)

        print "- Results written to: " + self.psv['project'] + ":" +run['resultsFolder']
        if self.multi_rep:
            for ltr in sorted( self.psv['reps'].keys() ):
                rep = self.psv['reps'][ltr]
                print "                      " + self.psv['project'] + ":" +rep['resultsFolder']

        print "- Steps to run:"
        to_run_count = 0
        if self.multi_rep:
            for ltr in sorted( self.psv['reps'].keys() ):
                rep = self.psv['reps'][ltr]
                for step in rep['path']:
                    if step in rep['stepsToDo']:
                        print "  * "+rep['rep_tech'] +': '+ rep['steps'][step]['app']+" will be run"
                    else:
                        if not step.find('concat') == 0:
                            print "    "+rep['rep_tech'] +': '+ rep['steps'][step]['app']+" has already been run"
                to_run_count += len(rep['stepsToDo'])
            if self.combined_reps:
                to_run_count += len(run['stepsToDo'])
                for step in run['path']:
                    if step in run['stepsToDo']:
                        print "  * "+run['rep_tech'] +': '+ run['steps'][step]['app']+" will be run"
                    else:
                        if not step.find('concat') == 0:
                            print "    "+run['rep_tech'] +': '+ run['steps'][step]['app']+" has already been run"
        else:
            to_run_count += len(run['stepsToDo'])
            for step in run['path']:
                if step in run['stepsToDo']:
                    print "  * "+run['steps'][step]['app']+" will be run"
                else:
                    if not step.find('concat') == 0:
                        print "    "+run['steps'][step]['app']+" has already been run"
        if to_run_count == 0:
            print "* All expected results are in the results folder, so there is nothing to do."
            print "  If this experiment/replicate needs to be rerun, then use the --force flag to "
            print "  rerun all steps; or remove suspect results from the folder before relaunching."
            sys.exit(0)

        if self.multi_rep:
            for ltr in sorted( self.psv['reps'].keys() ):
                rep = self.psv['reps'][ltr]
                if len(rep['deprecate']) > 0:
                    deprecated = rep['resultsFolder']+"/deprecated/"
                    if rep['resultsFolder'].endswith('/'):
                        deprecated = run['resultsFolder']+"deprecated/"
                    print "Will move "+str(len(rep['deprecate']))+" prior result file(s) to '" + deprecated+"'."
                    for fid in rep['deprecate']:
                        print "  " + dxencode.file_path_from_fid(fid)
        # single-rep and combined work both do run:
        if not self.multi_rep or self.combined_reps:
            if len(run['deprecate']) > 0:
                deprecated = run['resultsFolder']+"/deprecated/"
                if run['resultsFolder'].endswith('/'):
                    deprecated = run['resultsFolder']+"deprecated/"
                print "Will move "+str(len(run['deprecate']))+" prior result file(s) to '" + deprecated+"'."
                for fid in run['deprecate']:
                    print "  " + dxencode.file_path_from_fid(fid)


    def workflow_report_and_build(self,run,proj_id=None,template=False,test=True,verbose=False):
        '''Builds the multi-rep/combined-rep workflow from parts, reporting plans and returning wf.'''
        # NOT EXPECTED TO OVERRIDE
        # Report the plans
        if self.multi_rep:
            self.report_run_plans()
        else:
            self.report_run_plans(run)

        print "Checking for currently running analyses..."
        self.check_run_log(run['resultsFolder'], proj_id, verbose=True)

        # Move old files out of the way...
        if self.multi_rep:
            for rep in run['reps'].values():
                if len(rep['deprecate']) > 0 and not test:
                    deprecated = rep['resultsFolder']+"deprecated/"
                    print "Moving "+str(len(rep['deprecate']))+" "+rep['rep_tech']+" prior result file(s) to '"+ \
                                                                                        deprecated+"'..."
                    dxencode.move_files(rep['deprecate'],deprecated,proj_id)
        if not self.multi_rep or self.combined_reps:
            if len(run['deprecate']) > 0 and not test:
                deprecated = run['resultsFolder']+"deprecated/"
                print "Moving "+str(len(run['deprecate']))+" "+run['rep_tech']+" prior result file(s) to '"+ \
                                                                                    deprecated+"'..."
                dxencode.move_files(run['deprecate'],deprecated,proj_id)
        
        # Build the workflow...
        wf = None
        if self.multi_rep:
            for rep_id in sorted( run['reps'].keys() ):
                rep = run['reps'][rep_id]
                if len(rep['stepsToDo']) > 0: # Going to get noisy with multimple reps
                    if test:
                        print "Testing workflow assembly for "+rep['rep_tech']+"..."
                    else:
                        print "Assembling workflow for "+rep['rep_tech']+"..."
                    wf = self.create_or_extend_workflow(rep, rep_id, wf=wf, test=test,template=template)
        if not self.multi_rep or self.combined_reps:
            if len(run['stepsToDo']) > 0:
                if test:
                    print "Testing workflow assembly for "+run['rep_tech']+"..."
                else:
                    print "Assembling workflow for "+run['rep_tech']+"..."
                wf = self.create_or_extend_workflow(run, None, wf=wf, test=test,template=template)
            
        return wf


    def check_run_log(self,results_folder,proj_id,verbose=False):
        '''Checks for currently running jobs and will exit if found.'''
        # NOT EXPECTED TO OVERRIDE
        run_log_path = results_folder + '/' + dxencode.RUNS_LAUNCHED_FILE
        log_fids = dxencode.find_file(run_log_path,proj_id,multiple=True,recurse=False)
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


    def log_this_run(self,run_id,results_folder):
        '''Adds a runId to the runsLaunched file in resultsFolder.'''
        # NOT EXPECTED TO OVERRIDE
        # NOTE: DX manual lies?!  Append not possible?!  Then write new/delete old
        run_log_path = results_folder + '/' + dxencode.RUNS_LAUNCHED_FILE
        old_fids = dxencode.find_file(run_log_path,self.proj_id,multiple=True,recurse=False)
        new_fh = dxpy.new_dxfile('w',project=self.proj_id,folder=results_folder, \
                                                                  name=dxencode.RUNS_LAUNCHED_FILE)
        new_fh.write(run_id+' started:'+str(datetime.now())+'\n')
        if old_fids is not None:
            for old_fid in old_fids:
                with dxpy.open_dxfile(old_fid) as old_fh:
                    for old_run_id in old_fh:
                        new_fh.write(old_run_id+'\n')
            proj = dxpy.DXProject(self.proj_id)
            proj.remove_objects(old_fids)
            
        new_fh.close()

    def launch_pad(self,wf,run,ignition=False):
        '''Launches or just advertises preassembled workflow.'''
        # NOT EXPECTED TO OVERRIDE
        if wf == None:
            print "ERROR: failure to assemble workflow!"
            sys.exit(1)

        if ignition:
            print "Launch sequence initiating..."
            wf_run = wf.run({}, project=self.proj_id)
            if wf_run == None:
                print "ERROR: failure to lift off!"
                sys.exit(1)
            else:
                print "  We have liftoff!"
                wf_dict = wf_run.describe()
                self.log_this_run(wf_dict['id'],run['resultsFolder'])
                print "  Launched " + wf_dict['id']+" as '"+wf.name+"'"
        else:
            print "Workflow '" + wf.name + "' has been assembled in "+run['resultsFolder'] + \
                                                                        ". Manual launch required."

    def run(self):
        '''Runs launch from start to finish using command line arguments.'''
        # NOT EXPECTED TO OVERRIDE

        args = self.get_args()

        print "Retrieving pipeline specifics..."
        self.psv = self.pipeline_specific_vars(args)
        print "Running in project ["+self.proj_name+"]..."
        
        print "Building apps dictionary..."
        rep_steps, file_globs = self.assemble_steps_and_globs(self.REP_STEP_ORDER)
        combined_steps = None
        if self.combined_reps:
            combined_steps, file_globs = self.assemble_steps_and_globs(self.COMBINED_STEP_ORDER,True)

        # finding fastqs and prior results in a stadardized way
        self.finding_rep_inputs_and_priors(rep_steps,file_globs,args.test)

        #print "Checking for control files..."
        self.find_all_control_files(args.test)

        # Checking for combined-replicate inputs and priors...
        self.finding_combined_inputs_and_priors(combined_steps,file_globs,args.test)

        # finding pipeline specific reference files in a stadardized way
        self.find_all_ref_files()

        # deternine steps to run in a stadardized way
        self.determine_steps_needed(rep_steps,combined_steps,args.force)

        # Preperation is done. Now build up multi-rep workflow
        if self.multi_rep:
            run = self.psv
            wf = self.workflow_report_and_build(run,self.proj_id,test=args.test,template=(not args.run))

        else:            
            # Preperation is done. At this point on we either run rep 'a' or combined.
            if not self.combined_reps:
                run = self.psv['reps']['a']
            else:
                run = self.psv
            wf = self.workflow_report_and_build(run,self.proj_id,test=args.test,template=(not args.run))

        # Exit if test only
        if args.test:
            print "TEST ONLY - exiting."
            sys.exit(0)

        # Roll out to pad and possibly launch
        self.launch_pad(wf,run,ignition=args.run)
                
        print "(success)"

if __name__ == '__main__':
    '''Run from the command line.'''
    # EXPECTED to crash and burn on base class!
    launch = Launch()
    launch.run()

