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
#        - PIPELINE_BRANCHES[branch_id]["ORDER"] is the list of steps in the pipeline branch
#        - PIPELINE_BRANCHES[branch_id]["STEPS"] contain step definitions and enforces dependencies by:
#          'input' file tokens matching to 'result' file tokens of earlier steps.
#        - FILE_GLOBS is needed for locating result files from prior runs.
#      - Combined replicate processing is optional and depends upon PIPELINE_BRANCHES["COMBINED_REPS"], etc.
#        - The standard combination expected PIPELINE_BRANCH_ORDER[ "REP", "COMBINED_REPS" ].
#        - More complex combinations can be acheived by overriding add_combining_reps() in the derived class.
#        - Tying results from one branch to inputs into another relies on input file tokens ending with '_a', '_b' or '_set'.
#      - Control files are supported but most likely will require function overrides in derived
#        classes (e.g. rampageLaunch.py)
#      - In PIPELINE_BRANCHES[branch_id]["STEPS"], both 'inputs' and 'results' are key:value lists where the key is used
#        to attach results of one step to the inputs of another, while values must match dxapp.json names.
#        Further, any input key ending in '_set' is expected to be an array of files (e.g. "input": { "reads_set": "reads" }  
#        will expect one or more files to be found and used for dxapp.json input named "reads" ).
#      - The FILE_GLOBS are matched by input/result keys, not values.
#        (e.g. "input": { "pooled_bam": "bam" } matches file_glob "pooled_bam" and dxapp.json input named "bam".
#        This example would allow the same applet to be used on replicate bams and poold bams in the same pipeline.)
#
#   A design note about replicates in the launchers:
#      - "rep" objects should contain (almost) all information necesary to launch processing against them.
#        A few minor things common to the entire experiment are stored directly in the psv (pipeline specific variables) object
#      - Each "rep" object will get processed by one PIPELINE_BRANCH which determines which STEPS are run.
#      - All ENCODEd replicates ready for processing should get loaded into psv["reps"], keyed with letters: "a","b","c"...
#      - These "simple reps" have a "rep_tech" like "rep2_1" (for bio_rep=2, tech_rep=1).
#      - Combined replicates have a "rep_tech" like "reps1_1-2_1" (for combining simple reps rep1_1 and rep2_1).
#      - Combined replicates should be keyed (in sort order) after their "tributary" reps (e.g. "a","b","b-combines-a-b").
#      - The final combined replicate (aka "SEA" into which all rivers flow) is keyed with self.SEA_ID ("zzz")
#      - In the 'standard combination model' (which can be overridden in descendent classes): 
#        - Simple reps are processed by the "REP" pipeline branch.
#        - The one combined rep (aka "SEA") is processed by the "COMBINED_REPS" pipeline branch

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
    
    PIPELINE_BRANCH_ORDER = None # MUST replace in derived class. Example: [ "REP", "COMBINED_REPS" ]
    '''A pipeline is frequently made of branches that flow into each other, such as replicate level to combined replicate.'''
    
    PIPELINE_BRANCHES = None # MUST replace in derived class. 'ORDER' and "STEPS" are required for each branch,
    # Example:               # though simple "STEPS" may be discovered.
    #{    "REP": {
    #            "ORDER": [ "rep_app_one", ... ],
    #     # NOTE: In simple cases "STEPS" can be discovered from dnanexus apps (example: srnaLaunch.py)
    #     # HINT: Try to discover the "STEPS" usubg build_simple_steps(verbose) to print the json, then modify from there  
    #            "STEPS": {
    #                        "rep_app_one": {
    #                            "inputs":  { "reads":    "reads" },
    #                            "app":       "rep_app_1", 
    #                            "params":  { "nthreads": "nthreads"}, 
    #                            "results": { "rep_bam":  "bam_out" }
    #                        }, ...
    #    },
    #    "COMBINED_REPS": ...
    #}
    '''Each branch must define the 'steps' and their (artificially) linear order.'''

    # "STEPS" objects within PIPELINE_BRANCHES contain the following:
    #    "step-name":   { # Must match "ORDER" and may match dx applet
    #        "app":     "dx-applet-name", # Must match applet
    #        "params":  { "param_token": "dx_input_name", ...}, # Non-files: Token must match psv key.
    #        "inputs":  { "input_token": "dx_input_name", ...},  # Files: Token *may* match name
    #        "results": { "result_token": "dx_output_name", ...}} # Files: Token *may* match name
    # Note, tokens *may* be the same as dx names, but MUST MATCH other tokens to define dependencies as: 
    #     stepA:result_token == stepB:input_token
    # Also result_tokens must be in FILE_GLOB.keys()
    #
    # NOTE: Simple pipelines have no reused steps or dx name to param complications and the 
    #       "STEPS" objects can be generated directly from dx apps.  If unsure, print the json from
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
    ''' This is the default Genome that experiments are mapped to.'''
    
    STANDARD_LABELS = {
        'long-rna-seq':  { 'long': 'Long-RNA-seq',  'short': 'lrna',  'se_and_pe': True  }, 
        'small-rna-seq': { 'long': 'Small-RNA-seq', 'short': 'srna',  'se_and_pe': False }, 
        'rampage':       { 'long': 'Rampage',       'short': 'ramp',  'se_and_pe': False }, 
        'dnase-seq':     { 'long': 'DNase-seq',     'short': 'dnase', 'se_and_pe': True  },
    }
    '''Standard labelling requires exp_type specific labels.  This can be overridden in descendent classes.'''
    
    SEA_ID = 'zzz'
    '''SEA is the final branch into which all tributaries flow.'''

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

        if self.PIPELINE_BRANCH_ORDER != None and "COMBINED_REPS" in self.PIPELINE_BRANCH_ORDER:
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

    def add_combining_reps(self, psv):
        '''Defines how replicated are combined.'''
        # MAY NEED TO REPLACE for pipelines with combined steps
        
        reps = psv['reps']
        # In the 'standard combining model' PIPELINE_BRANCH_ORDER = [ "REP", "COMBINED_REPS" ]
        # and all ENCODEd replicates are in psv['reps'] keyed as 'a','b',... and having rep['rep_tech'] = 'rep1_1'
        # All these simple reps will have rep['branch_id'] = "REP"
        
        # If there are non-standard combinations, such as tech_reps combining, then bio_reps combining, define it here...
        # HINT: dnaseLaunch.py has a working example

        # If requested with -cr, OR if exactly 2 reps from different bio_reps were found then self.combined_reps == True
        # In the 'standard combining model' the one combined replicate keyed with self.SEA_ID is created below:
        # OVERRIDE as you will!
        if 'COMBINED_REPS' in self.PIPELINE_BRANCH_ORDER and self.combined_reps:
            assert len(reps.keys()) == 2
            assert 'a' in reps.keys()
            assert 'b' in reps.keys()
            
            sea = {} # SEA is the final branch into which all tributaries flow
            sea['branch_id'] = 'COMBINED_REPS'
            sea['tributaries'] = ['a','b']
            sea['rep_tech'] = 'reps' + reps['a']['rep_tech'][3:] + \
                                 '-' + reps['b']['rep_tech'][3:]  # combined delimited by '-'
            psv['reps'][self.SEA_ID] = sea
            psv['rep_tech'] = sea['rep_tech']  # psv gets labelled the same as the sea
            
        
    def find_combined_inputs(self,rep,steps,file_globs):
        '''Finds the inputs for a combined run in the directories of the replicate runs.'''
        # MAY NEED TO REPLACE for pipelines with combined steps
        
        if not self.combined_reps:
            return

        print "Checking for combined-replicate inputs..."
        inputs = {}
        for step in rep['path']:
            for file_token in steps[step]['inputs'].keys():
                if file_token[-2] != '_':
                    continue
                rep_key = file_token[-1].lower()
                if rep_key not in self.psv['reps'].keys():
                    print "*** " + rep_key + " not found in psv['reps']"
                    continue
                tributary = self.psv['reps'][rep_key]
                # TODO: No need for multiples at this time.  Deal with it when it comes up.
                fid = dxencode.find_file(tributary['resultsFolder'] + file_globs[file_token],\
                                                    self.proj_id, multiple=False, recurse=False)
                if fid != None:
                    rep['priors'][file_token] = fid
                    inputs[file_token] = [ fid ]
                elif not self.multi_rep:
                    print "Error: Necessary '%s' for combined run, not found in '%s'." \
                                                               % (file_token, rep['resultsFolder'])
                    print "       Please run for single replicate first."
                    sys.exit(1)
        
        if 'inputs' not in rep:
            rep['inputs'] = inputs
        else:
            rep['inputs'].update( inputs )
    

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
        
        self.exp = dxencode.get_exp(cv['experiment'],key=key)
        cv['exp_type'] = dxencode.get_assay_type(cv['experiment'],self.exp)
        
        # Load up the encoded and combining "reps" which will be processed by PIPELINE_BRANCHES
        self.load_reps(args, cv, cv['experiment'], self.exp)
        
        assert 'a' in cv['reps']

        # Only supported genomes
        cv['gender'] = cv['reps']['a']['sex']
        organism = cv['reps']['a']['organism']
        if organism in dxencode.GENOME_DEFAULTS:
            cv['genome'] = dxencode.GENOME_DEFAULTS[organism]
        else:
            print "Organism %s not currently supported" % organism
            sys.exit(1)

        # Paired ends?  Make sure combined reps reflect their tributaries and exp reflects all
        if self.combined_reps:
            cv['paired_end'] = False
            for rep_id in sorted( cv['reps'].keys() ):  # Must be in order
                if len(rep_id) == 1: # Simple reps are the tributaries and should have "paired_end" set already
                    assert "paired_end"      in cv['reps'][rep_id]
                    assert "tributaries" not in cv['reps'][rep_id]
                    continue
                river = cv['reps'][rep_id]
                assert "tributaries" in river
                paired = None
                for trib_id in river['tributaries']:
                    tributary = cv['reps'][trib_id]
                    if paired == None:
                        first_rep_tech = tributary['rep_tech']
                        paired = tributary['paired_end']
                        if 'paired_end' not in river:
                            river['paired_end'] = paired
                        if paired:
                            cv['paired_end'] = paired
                    elif paired != tributary['paired_end']:
                        print "- WARNING: Combining replicates '"+first_rep_tech+"' and '"+tributary['rep_tech']+"' " + \
                                        "are expected to be all paired-end or single-end!"
                        river['paired_end'] = False  # Mixes will be treated as single end!
                        #sys.exit(1)
        else:  # for independent simple replicates, top level doesn't really matter.
            cv['paired_end'] = cv['reps']['a']['paired_end']

        # Default locations
        cv['refLoc'] = args.refLoc
        if cv['refLoc'] == dxencode.REF_FOLDER_DEFAULT:
            cv['refLoc'] = dxencode.REF_FOLDER_DEFAULT + cv['genome'] + '/'
        if not cv['refLoc'].endswith('/'):
            cv['refLoc'] += '/' 
        cv['resultsLoc'] = dxencode.umbrella_folder(args.folder,self.FOLDER_DEFAULT,self.proj_name,cv['exp_type'],cv['genome'])
        cv['resultsFolder'] = cv['resultsLoc'] + cv['experiment'] + '/'
        self.update_rep_result_folders(cv)
        
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
                rep['name']  = cv['name']+"_"+rep['rep_tech']
            cv['title']   += " - "+cv['rep_tech']+" on "+cv['genome']+", "+cv['gender']
            cv['name']    += "_"+cv['rep_tech']

        return cv

    def load_reps(self,args, cv, exp_id, exp=None):
        '''
        Gathers replicates from encoded, then builds the cv['reps'] list from requested or all.
        '''
        # NOT EXPECTED TO OVERRIDE
        
        # A design note about replicates in the launchers:
        # - "rep" objects should contain (almost) all information necesary to launch processing against them.
        #   A few minor things common to the entire experiment are stored directly in the psv (pipeline specific variables)
        # - Each "rep" object will get processed by one PIPELINE_BRANCH which determines which STEPS are run.
        # - All ENCODEd replicates ready for processing should get loaded into psv["reps"], keyed with letters: "a","b","c"...
        # - These "simple reps" have a "rep_tech" like "rep2_1" (for bio_rep=2, tech_rep=1).
        # - Combined replicates have a "rep_tech" like "reps1_1-2_1" (for combining simple reps rep1_1 and rep2_1).
        # - Combined replicates should be keyed (in sort order) after their "tributary" reps (e.g. "a","b","b-combines-a-b").
        # - The final combined replicate (aka "SEA" into which all rivers flow) is keyed with self.SEA_ID ("zzz")
        # - In the 'standard combination model' (which can be overridden in descendent classes): 
        #   - Simple reps are processed by the "REP" pipeline branch.
        #   - The one combined rep (aka "SEA") is processed by the "COMBINED_REPS" pipeline branch
        cv_reps = {}
        full_mapping = dxencode.get_full_mapping(exp_id,exp=exp)
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
            if self.PIPELINE_BRANCH_ORDER != None and 'COMBINED_REPS' in self.PIPELINE_BRANCH_ORDER:
                if len(reps) == 2 and cv_reps['a']['br'] != cv_reps['b']['br']:
                    self.combined_reps = True
            self.multi_rep = (len(reps) > 1)
            
        # A little more rep tidying
        for ltr in cv_reps.keys():
            rep = cv_reps[ltr]
            rep['branch_id'] = 'REP'  # MUST override for non-standard pipelines which have ['REP','COMBINED_REPS']
            rep['concat_id'] = 'reads'
            if 'paired_end' in rep and rep['paired_end']: 
                rep['concat_id2'] = 'reads2'
                
        # mult-rep rep_tech: 
        if self.combined_reps:
            cv['rep_tech'] = 'reps' + cv_reps['a']['rep_tech'][3:] + \
                                '-' + cv_reps['b']['rep_tech'][3:]  # combined delimited by '-'
        elif self.multi_rep:
            cv['rep_tech'] = 'reps:' + cv_reps['a']['rep_tech'][3:]
            for ltr in sorted(cv_reps.keys()):
                if ltr != 'a':
                    cv['rep_tech'] += ','+cv_reps[ltr]['rep_tech'][3:] # not combined are delimited by ','
        else:
            cv['rep_tech'] = cv_reps['a']['rep_tech']
            
        cv['reps'] = cv_reps 

        # Call override-able function to define replicate combinations.
        self.add_combining_reps(cv)


    def update_rep_result_folders(self, cv):
        ''' create input object for a step and extends the file_globs dict as appropriate.'''
        # NOT EXPECTED TO OVERRIDE
        for rep_id in sorted( cv['reps'].keys() ):
            rep = cv['reps'][rep_id]
            if rep_id == self.SEA_ID:  # SEA is the ultimate branch where all pipelines flow
                rep['resultsFolder'] = cv['resultsFolder']
            elif 'tributaries' not in 'rep':
                rep['resultsFolder'] = cv['resultsFolder'] + rep['rep_tech'] + '/'
            elif 'resultsFolder' in 'rep':
                rep['resultsFolder'] = cv['resultsFolder'] + rep['resultsFolder'] + '/'
            else:
                rep['resultsFolder'] = cv['resultsFolder']


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


    def assemble_steps_and_globs(self,branch_id,branch):
        '''Assemble the requested steps and file globs.'''
        # NOT EXPECTED TO OVERRIDE
        assert 'ORDER' in branch
        pipe_path = branch['ORDER']    

        if self.FILE_GLOBS != None and 'STEPS' in branch:
            # Use for complex pipelines:
            file_globs = self.FILE_GLOBS
            pipe_steps = branch['STEPS']
        else:
            # Use for simple pipelines
            pipe_steps, file_globs = self.build_simple_steps(pipe_path,self.proj_id)

        for rep in self.psv['reps'].values():
            if rep["branch_id"] != branch_id:
                continue
            if type(pipe_path) == dict:
                assert 'se' in pipe_path and 'pe' in pipe_path
                if rep['paired_end']:
                    rep['path'] = pipe_path['pe']
                else:
                    rep['path'] = pipe_path['se']
            else:
                rep['path'] = pipe_path
            rep['steps'] = pipe_steps
            
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


    def find_inputs_and_priors(self,branch_id,globs,test):
        '''Finds the inputs and priors for all reps for a given branch of the pipeline.'''
        # NOT EXPECTED TO OVERRIDE
        
        print "Checking for prior results..."
        # NOTE: priors is a dictionary of fileIds that will be used to determine stepsToDo
        #       and fill in inputs to workflow steps
        for rep in self.psv['reps'].values():
            if rep['branch_id'] != branch_id:
                continue
            if not test:
                if not dxencode.project_has_folder(self.project, rep['resultsFolder']):
                    self.project.new_folder(rep['resultsFolder'],parents=True)
            rep['priors'] = self.find_prior_results(rep['path'],rep['steps'],rep['resultsFolder'],globs)

        print "Checking for input files..."
        # Find all reads files and move into place
        # TODO: files could be in: dx (usual), remote (url e.g.https://www.encodeproject.org/...)
        #       or possibly local, Currently only DX locations are supported.
        for rep_id in self.psv['reps'].keys():
            if len(rep_id) > 1:  # input fastqs only expected for simple ENCODEd reps (keyed by single letters)
                continue
            rep = self.psv['reps'][rep_id]
            if 'fastqs' not in rep:
                continue
            rep['inputs'] = {}
            reads_token = 'reads'
            if rep['paired_end']:
                reads_token += '1'
            rep['inputs']['Reads1'] = dxencode.find_and_copy_read_files(rep['priors'], \
                                                rep['fastqs']['1'],  test, reads_token, \
                                                rep['resultsFolder'], False, self.proj_id)
            # Note: rep['fastqs']['2'] and rep['inputs']['Reads2'] will be empty on single-end
            rep['inputs']['Reads2'] = dxencode.find_and_copy_read_files(rep['priors'], \
                                                rep['fastqs']['2'], test, 'reads2', \
                                                rep['resultsFolder'], False, self.proj_id)

        for rep in self.psv['reps'].values():
            # If rep has tributaries, then look for those priors!
            if rep['branch_id'] != branch_id:
                continue
            if 'tributaries' in rep:
                self.find_combined_inputs(rep,rep['steps'],globs)


    def find_all_ref_files(self):
        '''Locates all reference files based upon organism and gender.'''
        # NOT EXPECTED TO OVERRIDE
        
        print "Looking for reference files..."
        ref_priors = {}
        self.find_ref_files(ref_priors)
        
        # Add into all rep priors
        for rep in self.psv['reps'].values():
            rep['priors'].update(ref_priors)


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


    def determine_steps_needed(self, force=False):
        '''Determine steps needed for replicate(s) and combined, base upon prior results.'''
        # NOT EXPECTED TO OVERRIDE
        
        print "Determining steps to run..."
        # NOTE: stepsToDo is an ordered list of steps that need to be run
        for rep_id in sorted( self.psv['reps'].keys() ):
            rep = self.psv['reps'][rep_id]
            rep['deprecate'] = [] # old results will need to be moved/removed if step is rerun
            rep['stepsToDo'] = self.determine_steps_to_run(rep['path'], rep['steps'], \
                                                    rep['priors'], rep['deprecate'], force=force)


    def find_results_from_prev_branch(self, river, inp_token, expect_set):
        '''
        If combined_reps step, look for input from results of replicate level step.
        '''
        if not self.combined_reps:    # Not even doing a wf with combining branches 
            return None
        if 'tributaries' not in river:  # Not a combining branch
            return None
        if inp_token not in self.FILE_GLOBS: # assertable
            return None
        # The combined_step input tokens should either end in '_a', '_b' or the like,
        # Or the step expects to combine N files so expect_set=True and all matching results should be combined.
        
        # Using letters since 2 tech_rep branches could be combined into 2 bio_rep branches and then 1 experiment level branch
        # Thus tr branch_ids would be: 'a','b','c','d'
        #  and br branch_ids might be 'b-bio_rep1','d-bio_rep2', and those branches flow into the sea (self.SEA_ID='zzz').
        # So if a BRANCH STEP for combined rep 'b-bio_rep1' expects 2 inputs ending '_a' and '_b then
        #       the same STEP for combined rep 'd-bio_rep2' would expect the same 2 endings, not '_c' and '_d'
        letters = deque('abcdefghijklmnopqrstuvwxyz')
        looking_for = None
        if expect_set:
            # Expectation that when combining results, dx app's input names either end in '_a', etc. or '_set'
            if not inp_token.endswith('_set'):
                print "ERROR: set is expected but dx app's input name did not end in '_set'"
                return None
            results_array = []
        elif inp_token[-2] == '_' and inp_token.lower()[-1] in letters:
            looking_for = inp_token.lower()[-1]
        elif looking_for == None: 
            # Can we do this?  NOT an warning: no current means to distinguish between expecting a prior from not expecting one
            #print "- WARNING: Not looking in prior branch because DX app's input name does not end in '_a', '_b' or '_set' + \
            #                                                        " (inp_token:"+inp_token+" inp_token:"+inp_token+")."
            return None
            
        for tributary_id in river['tributaries']:
            ltr = letters.popleft()
            if looking_for != None and looking_for != ltr:
                continue
            assert tributary_id in self.psv['reps']
            tributary = self.psv['reps'][tributary_id]
            if 'prevStepResults' not in tributary:
                continue
            prev_results = tributary['prevStepResults']
            # How to match prevResults[inp_token='bwa_bam'] in prevResults to compare to step['inputs'][inp_token='bam_a']
            # The only way is by matching file globs which MUST tie branch leaps together.
            for result_token in prev_results.keys():
                if result_token in self.FILE_GLOBS:
                    # TODO: ultimately these are paired by matching GLOBS!  This frailty might be avoided
                    if self.FILE_GLOBS[result_token] == self.FILE_GLOBS[inp_token]:
                        if expect_set:
                            #print "- Expected set in branch priors and found inp:'"+inp_token + \
                            #                                                "' from prior branch result:'"+result_token+"'" 
                            results_array.append(prev_results[result_token])
                            break  # Done with this branch, but look for more results in next branches
                        if looking_for == ltr:
                            #print "- Expected '"+ltr+"' in branch priors and found inp:'"+inp_token + \
                            #                                                "' from prior branch result:'"+result_token+"'" 
                            return prev_results[result_token] # found result on branch matched by letter
                        #elif looking_for == None:
                        #    return prev_results[result_token] # return first match from first branch
        
        # Went through all branches so return any results found
        if not expect_set or len(results_array) == 0:
            return None
        return results_array


    def create_or_extend_workflow(self,rep, rep_id, wf=None,app_proj_id=None,test=False,template=False):
        '''
        This function will populate or extend a workflow with the steps in rep['stepsToDo'] and return 
        the workflow unlaunched.   It relies on the rep object containing *almost* all information to run
        the PIPELINE_BRANCH associated with the rep.  The rep should contain rep['priors'] defining existing files;
        rep['steps'] objects defining applet, input and output requirements; etc.  A few common parameters will
        be pulled from the psv (pipeline specific variables) object.
        '''
        # NOT EXPECTED TO OVERRIDE

        if len(rep['stepsToDo']) < 1:
            return None
        if app_proj_id == None:
            app_proj_id = self.proj_id

        # create a workflow object
        if self.multi_rep:
            wf_folder = self.psv['resultsFolder']
            wf_name = self.psv['name']
        else:
            wf_name = rep['name']
            wf_folder = rep['resultsFolder']
        if not test and wf == None:
            wf = dxpy.new_dxworkflow(title=wf_name,name=wf_name,folder=wf_folder,
                                           project=self.proj_id,description=self.psv['description'])

        # NOTE: prevStepResults dict contains links to result files to be generated by previous steps
        if 'prevStepResults' not in rep:
            rep['prevStepResults'] = {}
        steps = rep['steps']
        for step in rep['stepsToDo']:
            appName = steps[step]['app']
            app = dxencode.find_applet_by_name(appName, app_proj_id)
            inp_defs = app.describe().get('inputSpec') or []
            appInputs = {}
            # file inputs
            for fileToken in steps[step]['inputs'].keys():
                appInp = steps[step]['inputs'][fileToken]
                # Need inp_def to see if set is expected
                for inp_def in inp_defs:
                    if inp_def["name"] == appInp:
                        break
                if inp_def["name"] != appInp:
                    print "ERROR: Pipeline definition for applet '"+appName+"' fileToken '"+fileToken+"'"
                    print "value of '"+appInp+"' not found in DX definition."
                    sys.exit(1)
                expect_set = (inp_def["class"] == 'array:file')
                
                # Now try to find prior results or existing files to fill in this input
                prev_results = None
                if 'tributaries' in rep:
                    # First: check in previous branch (e.g. rep level pipeline) is input to the combining step
                    prev_results = self.find_results_from_prev_branch(rep, fileToken, expect_set)
                if prev_results != None:
                    appInputs[ appInp ] = prev_results
                elif fileToken in rep['prevStepResults']:
                    # Now the most likely: previous results of this branch
                    appInputs[ appInp ] = rep['prevStepResults'][fileToken]
                elif fileToken in rep['priors'] \
                  or (expect_set and fileToken + "_set" in rep['priors']):
                    # Finally, look in prior files found in results folder
                    if expect_set:
                        alt_token = fileToken
                        if fileToken not in rep['priors']:
                            alt_token = fileToken + "_set"  # Sometimes sets expect  the token to end in '_set'!
                        if isinstance(rep['priors'][alt_token], list):
                            #print "- Expecting input set and found "+str(len(rep['priors'][alt_token]))+" file(s)."
                            appInputs[ appInp ] = []
                            for fid in rep['priors'][alt_token]:
                                appInputs[ appInp ] += [ dxencode.FILES[fid] ]
                        else:
                            appInputs[ appInp ] = [ dxencode.FILES[ rep['priors'][fileToken] ] ]
                    else:
                        assert(not isinstance(rep['priors'][fileToken], list))
                        appInputs[ appInp ] = dxencode.FILES[ rep['priors'][fileToken] ]
                #elif not template:
                #    print "- Looking for existing files but '"+fileToken+"' not in rep['priors']!"
                if appInp not in appInputs and not template:
                    print "ERROR: step '"+step+"' can't find input '"+fileToken+"'!"
                    sys.exit(1)
            
            # Non-file app inputs
            if 'params' in steps[step]:
                for param in steps[step]['params'].keys():
                    appParam = steps[step]['params'][param]
                    if param in rep:
                        appInputs[ appParam ] = rep[param]
                    elif param in self.psv:
                        appInputs[ appParam ] = self.psv[param]
                    elif not template:
                        print "ERROR: step '"+step+"' unable to locate '"+param+ \
                                                          "' in pipeline specific variables (psv)."
                        sys.exit(1)
            
            # Now we are ready to add wf stage
            if not test:
                stageId = wf.add_stage(app, stage_input=appInputs, folder=rep['resultsFolder'])
            # outputs, which we will need to link to
            for fileToken in steps[step]['results'].keys():
                appOut = steps[step]['results'][fileToken]
                if test:
                    rep['prevStepResults'][ fileToken ] = 'fake-for-testing'
                else:
                    rep['prevStepResults'][ fileToken ] = dxpy.dxlink({'stage': stageId, \
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
        print "- Inputs:"
        if self.multi_rep:
            for ltr in sorted( self.psv['reps'].keys() ):
                rep = self.psv['reps'][ltr]
                rep_tech_msg = rep['rep_tech']
                for input_type in sorted( rep['inputs'].keys() ):
                    input_type_msg = input_type
                    if len(rep['inputs'][input_type]) > 0: 
                        for fid in rep['inputs'][input_type]:
                            print "  - "+rep_tech_msg+' '+input_type_msg+' '+dxencode.file_path_from_fid(fid)
                            rep_tech_msg = '      '  # Trick to make output less cluttered
                            input_type_msg = '      '
        else:
            rep_tech_msg = run['rep_tech']
            for input_type in sorted( run['inputs'].keys() ):
                input_type_msg = input_type
                if len(run['inputs'][input_type]) > 0: 
                    for fid in run['inputs'][input_type]:
                        print "  - "+rep_tech_msg+' '+input_type_msg+' '+dxencode.file_path_from_fid(fid)
                        rep_tech_msg = '      '
                        input_type_msg = '      '

        print "- Reference files:"
        # Should be enough to show the ref priors from 'a' rep for single-rep, multi or combined
        for token in self.psv['ref_files']:
            print "  " + dxencode.file_path_from_fid(self.psv['reps']['a']['priors'][token],True)
        
        print "- Results written to: " + self.psv['project'] + ":" +run['resultsFolder']
        if self.multi_rep:
            for rep_id in sorted( self.psv['reps'].keys() ):
                if rep_id != self.SEA_ID:
                    rep = self.psv['reps'][rep_id]
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
        else:
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
            for rep_id in sorted( self.psv['reps'].keys() ):
                rep = self.psv['reps'][rep_id]
                if len(rep['deprecate']) > 0 and not test:
                    deprecated = rep['resultsFolder']+"deprecated/"
                    print "Moving "+str(len(rep['deprecate']))+" "+rep['rep_tech']+" prior result file(s) to '"+ \
                                                                                        deprecated+"'..."
                    dxencode.move_files(rep['deprecate'],deprecated,proj_id)
        if not self.multi_rep:
            if len(run['deprecate']) > 0 and not test:
                deprecated = run['resultsFolder']+"deprecated/"
                print "Moving "+str(len(run['deprecate']))+" "+run['rep_tech']+" prior result file(s) to '"+ \
                                                                                    deprecated+"'..."
                dxencode.move_files(run['deprecate'],deprecated,proj_id)
        
        # Build the workflow...
        wf = None
        if self.multi_rep:
            for rep_id in sorted( self.psv['reps'].keys() ):
                rep = self.psv['reps'][rep_id]
                if len(rep['stepsToDo']) > 0: # Going to get noisy with multiple reps
                    if test:
                        print "Testing workflow assembly for "+rep['rep_tech']+"..."
                    else:
                        print "Assembling workflow for "+rep['rep_tech']+"..."
                    wf = self.create_or_extend_workflow(rep, rep_id, wf=wf, test=test,template=template)
        if not self.multi_rep:
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
        
        print "Building apps dictionaries..."
        # The point of using 'branches' instead of 'self.PIPELINE_BRANCHES' is that the app dicts *may* be discovered from dx
        branches = {}
        file_globs = {}
        for branch_id in self.PIPELINE_BRANCH_ORDER:
            steps_and_globs = self.assemble_steps_and_globs(branch_id,self.PIPELINE_BRANCHES[branch_id])
            branch = {}
            branch["steps"] = steps_and_globs[0]
            file_globs.update( steps_and_globs[1] )  # Should only need one set of file globs for all pipeline branches/steps
            branches[branch_id] = branch

        # finding fastqs, inputs and prior results in a stadardized way
        for branch_id in self.PIPELINE_BRANCH_ORDER:
            self.find_inputs_and_priors(branch_id,file_globs,args.test)

        # finding pipeline specific reference files in a stadardized way
        self.find_all_ref_files()

        # deternine steps to run in a stadardized way
        self.determine_steps_needed(args.force)

        # Preperation is done. Now build up multi-rep workflow
        if self.multi_rep or self.combined_reps:
            run = self.psv
        else:            
            run = self.psv['reps']['a']
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

