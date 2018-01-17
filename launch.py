#!/usr/bin/env python2.7
# launch.py 1.0.1

import argparse, os, sys, subprocess, json
from datetime import datetime
from collections import deque

import dxpy
#import dxencode
import dx
import encd

### TODO:
# 1) NEED TO MAKE a --template version not relying on ENCODEd at all!
# - Nice to have option to dx build pipeline applets
# - Nice to have Template module that mirrors Launch but resides in pipeline repo and does not use encd.py
# - Template should override exp, not require read files, not require ref files, not make folders.
# - Options: SE, PE, single/double replicate.  Example:
# - ./lrnaLaunch.py --template --rep 1 --se --use_refs --build_apps

# NOTES: This command-line utility will run the a pipeline for a single replicate or combined reps.
#      - All results will be written to a folder /<results_folder>/<exp_id>/rep<#>_<#>.
#      - If any results are already in that directory, then the steps that created those results
#        will not be rerun.
#      - If any jobs for the experiment and replicate are already running, nothing new will be
#        launched.
#      - Almost all the code is generic and derived classes will only need to extend or override a
#        a few things.  It relies upon JSON pipeline descriptions which in simple cases can be
#        directly read from dx (e.g. srnaLaunch.py).  More complex examples (e.g. lrnaLaunch.py)
#        must hard-code the pipeline json to handle repeated apps and dx name conflicts.
#        Tokens are used to abstract dx.app input/outout file names to avoid collisions.
#        - PIPELINE_BRANCHES[branch_id]["ORDER"] is the list of steps in the pipeline branch
#        - PIPELINE_BRANCHES[branch_id]["STEPS"] contain step definitions and enforces dependencies by:
#          'input' file tokens matching to 'result' file tokens of earlier steps.
#        - FILE_GLOBS is needed for locating result files from prior runs.
#      - Combined replicate processing is optional and depends upon PIPELINE_BRANCHES["COMBINED_REPS"], etc.
#        - The standard combination expected PIPELINE_BRANCH_ORDER[ "REP", "COMBINED_REPS" ].
#        - More complex combinations can be acheived by overriding add_combining_reps() in the derived class.
#        - Tying results from one branch to inputs into another relies on input file tokens ending with '_A', '_B' or '_ABC'.
#      - Control files are supported but most likely will require function overrides in derived
#        classes (e.g. rampageLaunch.py)
#      - In PIPELINE_BRANCHES[branch_id]["STEPS"], both 'inputs' and 'results' are key:value lists where the key is used
#        to attach results of one step to the inputs of another, while values must match dxapp.json names.
#      - The FILE_GLOBS are matched by input/result keys, not values.
#        (e.g. "input": { "pooled_bam": "bam" } matches file_glob "pooled_bam" and dxapp.json input named "bam".
#        This example would allow the same applet to be used on "replicate bams" and "poold bams" in the same pipeline.)
#      - RARE cases using 'output_values' and 'param_links' can support non-file outputs as inputs to later steps,
#        even across rep branch boundaries.
#
#   A design note about replicates in the launchers:
#      - "rep" objects should contain (almost) all information necesary to launch processing against them.
#        A few minor things common to the entire experiment are stored directly in the psv (pipeline specific variables) object
#      - Each "rep" object will get processed by one PIPELINE_BRANCH which determines which STEPS are run.
#      - All ENCODEd replicates ready for processing should get loaded into psv["reps"], keyed with letters: "a","b","c"...
#      - These "simple reps" have a "rep_tech" like "rep2_1" (for bio_rep=2, tech_rep=1).
#      - Combined replicates have a "rep_tech" like "reps1_1-1_2" (for combining simple reps rep1_1 and rep1_2).
#      - Combined replicates should be keyed (in sort order) after their "tributary" reps (e.g. "a","b","b-combines-a-b").
#      - The final combined replicate (aka "SEA" into which all rivers flow) is keyed with self.SEA_ID ("zzz")
#      - In the 'standard combination model' (which can be overridden in descendent classes):
#        - Simple reps are processed by the "REP" pipeline branch.
#        - The one combined rep (aka "SEA") is processed by the "COMBINED_REPS" pipeline branch
#
# Example derived classes:
#    https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/dnanexus/lrnaLaunch.py
#    https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/dnanexus/small-rna/srnaLaunch.py
#    https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/dnanexus/rampage/rampageLaunch.py
#    https://github.com/ENCODE-DCC/dna-me-pipeline/blob/master/dnanexus/dmeLaunch.py
#    https://github.com/ENCODE-DCC/dnase_pipeline/blob/master/dnanexus/dnaseLaunch.py

class Launch(object):
    '''
    Launch module is the common code for pipeline specific launcher classes.
    Launchers dynamically create and launch a workflow covering any needed steps of a pipeline.
    Steps needed are thos that have available inputs but not already generatedresults.
    '''

    PIPELINE_NAME = "MUST REPLACE IN DREIVED CLASS"
    ''' This must match the assay type returned by encd.get_assay_type({exp_id}).'''

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
    #        "results": { "result_token": "dx_output_name", ...}, # Files: Token *may* match name
    #        "output_values": { "output_token": "dx_output_name" }, # RARE: non-file output used as input to another step
    #        "param_links": { "param_token": {"name":"output_token", "rep":"sister"} }} # RARE: link input to non-file output
    # Note, tokens *may* be the same as dx names, but MUST MATCH other tokens to define dependencies as:
    #     stepA:result_token == stepB:input_token
    # Also result_tokens must be in FILE_GLOB.keys()
    #
    # NOTE: Simple pipelines have no reused steps or dx name to param complications and the
    #       "STEPS" objects can be generated directly from dx apps.  If unsure, print the json from
    #       self.build_simple_steps(verbose=True) to determine if it will work for your pipeline.

    PRUNE_STEPS = None
    '''Rare setting that allows dynamically pruning steps from a pipeline befpre running it.  Used by lrnaLaunch.py.'''

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
        'rampage':       { 'long': 'Rampage',       'short': 'ramp',  'se_and_pe': True },
        'dnase-seq':     { 'long': 'DNase-seq',     'short': 'dnase', 'se_and_pe': True  },
        'dna-me':        { 'long': 'DNA-methlyation','short': 'dme',  'se_and_pe': True  },
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
        self.template = False
        self.build_apps = False
        self.no_refs = False
        self.server_key = self.SERVER_DEFAULT # TODO: replace with self.encd.server_key when Encd class is created
        self.proj_name = None
        self.project = None
        self.proj_id = None
        self.exp = {}  # Will hold the encoded exp json
        self.psv = {} # will hold pipeline specific variables.
        self.multi_rep = False       # This run includes more than one replicate (e.g. rep1_1 and (rep2_1 or rep1_2)
        self.combined_reps = False   # This run includes 2 or more reps that flow into one
        self.combine_one_or_more = False  # Special case where one or more flow into one
        self.compare_techreps = False# Only for special cases where there is only 1 bio_rep but 2 techreps.
        self.detect_umi = False      # Only in DNase is there a umi setting buried in the fastq metadata.
        self.link_later = None       # Rare: when 2 sister branches link to each other, it requires a final wf pass.
        self.deprecate = []          # Master list of files to deprecate.  Needed to cross rep boundaries in steps_to_run
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
                        required=False)

        ap.add_argument('-r','--reps','--replicates',
                        help="Request one or more replicates (e.g.'rep1_1 2 2_2').'",
                        nargs='+',
                        required=False)

        ap.add_argument('-br','--br', '--biological-replicate',
                        help="Request single biological replicate number (--tr==1 implied)",
                        type=int,
                        #default='1',
                        required=False)

        ap.add_argument('-tr','--tr', '--technical-replicate',
                        help="Select single technical replicate number (default: 1, must select --br as well)",
                        type=int,
                        default='1',
                        required=False)

        if self.CONTROL_FILE_GLOB != None:
            # Include this argument for pipelines with controls
            ap.add_argument('-c', '--control',
                            help='The control file for this experiment (may comma delimit rep1,rep2).',
                            default=None,
                            required=False)

            ap.add_argument('-cp', '--control_path',
                            help="Path to look in for control files (default: '" + self.CONTROL_ROOT_FOLDER + "')",
                            default=self.CONTROL_ROOT_FOLDER,
                            required=False)

        ap.add_argument('--project',
                        help="Project to run analysis in (default: '" + \
                                                        dx.env_get_current_project() + "')",
                        required=False)

        ap.add_argument('--refLoc',
                        help="The location to find reference files (default: '" + \
                            self.REF_PROJECT_DEFAULT + ":" + self.REF_FOLDER_DEFAULT + "')",
                        default=self.REF_FOLDER_DEFAULT,
                        required=False)

        ap.add_argument('-f','--folder',
                        help="The location to to place results folders (default: '<project>:" + \
                                                                  self.FOLDER_DEFAULT + "')",
                        default=self.FOLDER_DEFAULT,
                        required=False)

        ap.add_argument('-g','--genome',
                        help="The genome assembly to run on (default: '<project>:" + \
                                                                  self.GENOME_DEFAULTS['human'] + "')",
                        default=None,
                        required=False)

        ap.add_argument('-ct', '--compare_techreps',
                        help='Only for special cases where there is only 1 bio_rep but 2 techreps.',
                        action='store_true',
                        required=False)

        ap.add_argument('--run',
                        help='Run the workflow after assembling it.',
                        action='store_true',
                        required=False)

        ap.add_argument('--server',
                        help="Server to post files to (default: '" + self.SERVER_DEFAULT + "')",
                        default=self.SERVER_DEFAULT,
                        required=False)

        ap.add_argument('--build_apps',
                        help='Build applets using build_apps script before building workflow.',
                        action='store_true',
                        required=False)

        ap.add_argument('--template',
                        help='Build a template workflow only.',
                        action='store_true',
                        required=False)

        ap.add_argument('--pe',
                        help='For template only: make paired-end mapping template.',
                        action='store_true',
                        required=False)

        ap.add_argument('--no_refs',
                        help='For template only: do not use standard ENCODE reference files as inputs.',
                        action='store_true',
                        required=False)

        ap.add_argument('--test',
                        help='Test run only, do not launch anything.',
                        action='store_true',
                        required=False)

        ap.add_argument('--force',
                        help='Force rerunning all steps.',
                        action='store_true',
                        required=False)

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
        psv = self.common_variables(args)

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
            filePath = self.psv['refLoc']+self.REFERENCES_FILES[ref][self.psv['genome']][self.psv['gender']]
            fid = self.find_file(filePath,self.REF_PROJECT_DEFAULT)
            if fid == None:
                sys.exit("ERROR: Unable to locate ref file: '" + filePath + "'")
            else:
                priors[ref] = fid

        # Or explicitly by key
        chrom_sizes = self.psv['refLoc']+self.REFERENCE_FILES['chrom_sizes'][self.psv['genome']][self.psv['gender']]
        chrom_sizes_fid = self.find_file(chrom_sizes,self.REF_PROJECT_DEFAULT)
        if chrom_sizes_fid == None:
            sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chrom_sizes + "'")
        else:
            priors['chrom_sizes'] = chrom_sizes_fid
        self.psv['ref_files'] = self.REFERENCE_FILES.keys()
        return priors


    def find_control_file(self,rep,default=None):
        '''Attempts to find an appropriate control file.'''
        # MAY NEED TO REPLACE for pipelines with controls

        if self.CONTROL_FILE_GLOB == None or self.template:
            return None
        if 'controls' not in rep:
            return None
        for (control_exp_id,control_rep_tech) in rep['controls']:
            control_root = self.psv['control_path']
            # FIXME: Cheating:
            if self.proj_name == "scratchPad" and self.psv['control_path'] == self.CONTROL_ROOT_FOLDER:
                control_root = "/lrna"
            path_n_glob = control_root + control_exp_id + '/' + control_rep_tech + '/' + self.CONTROL_FILE_GLOB
            target_folder = dx.find_folder(control_exp_id + '/' + control_rep_tech,self.project,control_root)
            #print >> sys.stderr, "=== Target folder [%s]" % target_folder
            if target_folder != None:
                path_n_glob = target_folder + self.CONTROL_FILE_GLOB
                #print >> sys.stderr, "=== path_n_glob [%s] Project:[%s]" % (path_n_glob,self.proj_id)
            fid = self.find_file(path_n_glob,self.proj_id,multiple=False,recurse=False)
            if fid != None:
                return dx.file_path_from_fid(fid)

        if default != None:
            return default
        #print json.dumps(rep,indent=4)
        print >> sys.stderr, "ERROR: Unable to find control in search of %s" % rep['controls']
        sys.exit(1)


    def find_all_control_files(self):
        '''Attempts to find all control files needed for a run
        .'''
        # MAY NEED TO REPLACE for pipelines with controls

        if self.CONTROL_FILE_GLOB == None or self.template:
            return
        print "Checking for control files..."
        def_controls = []
        if 'control' in self.psv and isinstance(self.psv['control'],str):
            def_controls = self.psv['control'].split(',')
        ix = 0
        for rep in self.psv['reps'].values():
            control = None
            if len(def_controls) > ix:
                control = self.find_control_file(rep,def_controls[ix])
            else:
                control = self.find_control_file(rep,None)
            ix += 1
            if control != None:
                rep['inputs']['Control'] = dx.find_and_copy_read_files(rep['priors'], \
                                                    [ control ], self.test, 'control_bam', \
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
            # Special case of 2 allows for designating sisters
            psv['reps'][sea['tributaries'][0]]['sister'] = sea['tributaries'][1]
            psv['reps'][sea['tributaries'][1]]['sister'] = sea['tributaries'][0]
            psv['reps'][self.SEA_ID] = sea
            psv['rep_tech'] = sea['rep_tech']  # psv gets labelled the same as the sea


    def find_combined_inputs(self,rep,steps,file_globs,verbose=False):
        '''Finds the inputs for a combined run in the directories of the replicate runs.'''
        # MAY NEED TO REPLACE for pipelines with combined steps
        #verbose=True

        if not self.combined_reps or 'tributaries' not in rep:
            return

        print "Detecting combined-replicate inputs for %s..." % rep['rep_tech']
        inputs = {}
        for step in rep['path']:
            if verbose:
                print >> sys.stderr, "DEBUG: step: "+ step
            for file_token in steps[step]['inputs'].keys():
                if verbose:
                    print >> sys.stderr, "DEBUG:   file_token: "+ file_token
                if file_token[-2] == '_':
                    letter = file_token[-1].lower()
                    if letter in rep['tributaries']:
                        tributaries = [ letter ]
                    elif rep['tributaries'][0][1] == '-': # second level combining? (e.g. tributaties = ['b-biorep','e-biorep'])
                        # Difficulty is that inp_a/inp_b are being matched to ['b-biorep','e-biorep'] by position
                        letters = 'abcdefghijklmnopqrstuvwxyz'
                        ix = letters.find(letter) # hopefully returns 0, then 1
                        if ix != -1 and ix < len(rep['tributaries']):
                            tributaries = [ rep['tributaries'][ix] ]  # Too clever by half?
                    if len(tributaries) == 0:
                        print "*** " + rep_key + " not found in rep['tributaries']"
                        continue
                elif file_token.endswith('_ABC'):  # Also too clever by 3/4ths: ending says all tributaries contribute to set.
                    tributaries = rep['tributaries']
                else:
                    continue
                if verbose:
                    print >> sys.stderr, "DEBUG:   - found "+str(len(tributaries))+" tributaries."
                rep['priors'][file_token] = []
                inputs[file_token] = []
                if not self.template:
                    for rep_key in tributaries:
                        if verbose:
                            print >> sys.stderr, "DEBUG:   - tributary: " + rep_key
                        tributary = self.psv['reps'][rep_key]
                        fid = self.find_file(tributary['resultsFolder'] + file_globs[file_token],\
                                                        self.proj_id, multiple=False, recurse=False)
                        if fid != None:
                            if len(tributaries) == 1:
                                rep['priors'][file_token] = fid
                            else:
                                rep['priors'][file_token].append( fid )
                            inputs[file_token].append( fid )
                        elif not self.multi_rep and not self.combine_one_or_more:
                            print >> sys.stderr, "ERROR: Necessary '%s' for combined run, not found in '%s'." \
                                                                       % (file_token, rep['resultsFolder'])
                            print >> sys.stderr, "       Please run for single replicate first."
                            sys.exit(1)

        if inputs:
            if verbose:
                print >> sys.stderr, "DEBUG: inputs:"
                print >> sys.stderr, json.dumps(inputs,indent=4,sort_keys=True)
                #sys.exit(1)
        if 'inputs' not in rep:
            rep['inputs'] = inputs
        else:
            rep['inputs'].update( inputs )


    ############## WRAPPED DXENCODE METHODS (makes templating without encd/dx easier) ############
    # Wrap a few global defines to make inheritance easier
    REF_PROJECT_DEFAULT = dx.REF_PROJECT_DEFAULT
    REF_FOLDER_DEFAULT = dx.REF_FOLDER_DEFAULT
    GENOME_DEFAULTS = dx.GENOME_DEFAULTS

    def find_file(self,filePath,project=None,verbose=False,multiple=False, recurse=True):
        '''Using a DX style file path, find the file.'''
        # wrap the dx version to make templating easier
        return dx.find_file(filePath,project=project,verbose=verbose,multiple=multiple,recurse=recurse)

    def umbrella_folder(self,folder,default,proj_name=None,exp_type=None,genome=None,annotation=None):
        '''Returns a normalized umbrella folder (that holds the experiments of a given type).'''
        if self.no_refs: # (no_refs is only True when templating)
            genome = None # If templating with no refs then this will hide genome and annotation
        # NOTE: annotation is no longer used in organizing folders.
        return dx.umbrella_folder(folder,default,proj_name,exp_type=exp_type,sub_folder="runs/",genome=genome)

    ############## NOT EXPECTED TO OVERIDE THE FOLLOWING METHODS ############
    def common_variables(self,args,fastqs=True):
        '''Initializes dict with common variables from args and encoded.'''
        # NOT EXPECTED TO OVERRIDE

        self.server_key = args.server
        encd.set_server_key(self.server_key) # TODO: change to self.encd = Encd(self.server_key)
        self.test = args.test
        if args.experiment == None and not args.template:
            print >> sys.stderr, "ERROR: either --experiment or --template is required."
            sys.exit(1)
        if args.template:
            self.template = True
            if args.no_refs:
                self.no_refs = True
        if args.build_apps:
            self.build_apps = True
        if args.compare_techreps:
            self.compare_techreps = True

        cv = {}
        self.proj_name = dx.env_get_current_project()
        if self.proj_name == None or args.project != None:
            self.proj_name = args.project
        if self.proj_name == None:
            print >> sys.stderr, "ERROR: Please enter a '--project' to run in."
            sys.exit(1)
        self.project = dx.get_project(self.proj_name)
        self.proj_id = self.project.get_id()

        cv['project']    = self.proj_name
        cv['experiment'] = args.experiment
        if args.experiment == None and self.template:
            cv['experiment'] = 'template'

        # A few experiment types have controls
        if self.CONTROL_FILE_GLOB != None:
            if args.control != None: # Very few cases do you want to explicitly set this
                cv['control'] = args.control
            cv['control_path'] = args.control_path  # My defining this, you can shorten the search time through dx paths

        if not self.template:
            print "Retrieving experiment specifics..."
            self.exp = encd.get_exp(cv['experiment'])
            cv['exp_type'] = encd.get_assay_type(cv['experiment'],self.exp)
            if cv['exp_type'] != self.PIPELINE_NAME:
                print >> sys.stderr, "ERROR: Experiment %s is not for '%s' but for '%s'" \
                                               % (cv['experiment'],self.PIPELINE_NAME,cv['exp_type'])
                sys.exit(1)
            if 'internal_status' in self.exp and self.exp['internal_status'] in encd.INTERNAL_STATUS_BLOCKS:
                print "ERROR: Experiment %s with internal_status of '%s' cannot be launched." % \
                                                                        (cv['experiment'],self.exp['internal_status'])
                sys.exit(1)

            lab = self.exp.get('lab','')
            if isinstance(lab,dict):
                cv["lab"] = lab.get('name','unidentified')
            else:
                cv["lab"] = lab[6:-1] # e.g. "/labs/gregory-crawford/"
            self.load_reps(args, cv, cv['experiment'], self.exp)
        else:
            print "Templating experiment specifics..."
            cv['exp_type'] = self.PIPELINE_NAME  # By definition, this should be correct
            self.load_template_reps(args, cv)

        assert 'a' in cv['reps']

        # Only supported genomes
        cv['gender'] = cv['reps']['a']['sex']
        organism = cv['reps']['a']['organism']
        if organism in self.GENOME_DEFAULTS:
            cv['genome'] = self.GENOME_DEFAULTS[organism]
            if args.genome != None:
                if args.genome != self.GENOME_DEFAULTS['human']:
                    if args.genome not in self.GENOMES_SUPPORTED:
                        print >> sys.stderr, "ERROR: Genome %s not currently supported" % args.genome
                        sys.exit(1)
                elif cv['genome'] != args.genome:
                    print >> sys.stderr, "ERROR: Genome %s does not match discovered genome %s" % (args.genome, cv['genome'])
                    sys.exit(1)
                cv['genome'] = args.genome
        else:
            print >> sys.stderr, "ERROR: Organism %s not currently supported" % organism
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
                        if cv['exp_type'] == 'dnase-seq':
                            print >> sys.stderr, "- Combining replicates '"+first_rep_tech+"' and '"+tributary['rep_tech']+"' " + \
                                            "are paired-end and single-end.  Demoting bio_rep to SE!"
                            river['paired_end'] = False  # Mixes will be treated as single end!
                        else:
                            print >> sys.stderr, "- ERROR: Combining replicates '"+first_rep_tech+"' and '"+tributary['rep_tech']+"' " + \
                                            "are expected to be all paired-end or single-end!"
                            sys.exit(1)
        else:  # for independent simple replicates, top level doesn't really matter.
            cv['paired_end'] = cv['reps']['a']['paired_end']

        # Special case for RNA pipelines:
        if cv['exp_type'] in ["long-rna-seq", "small-rna-seq", "rampage", "cage"]:
            (cv['stranded'], cv['strand_direction']) = encd.is_stranded(self.exp)
        # Special case for lrna paired_end:
        if cv['paired_end'] and cv['exp_type'] == "long-rna-seq":
            cv['ScriptSeq'] = encd.is_script_seq(self.exp)

        # Default locations
        cv['refLoc'] = args.refLoc
        if cv['refLoc'] == self.REF_FOLDER_DEFAULT:
            cv['refLoc'] = self.REF_FOLDER_DEFAULT
            if not self.no_refs:
                cv['refLoc'] += cv['genome'] + '/'
        elif not cv['refLoc'].endswith('/'):
            cv['refLoc'] += '/'
        genome = cv['genome']
        cv['resultsLoc'] = self.umbrella_folder(args.folder,self.FOLDER_DEFAULT,self.proj_name,cv['exp_type'],genome)
        cv['resultsFolder'] = cv['resultsLoc']
        if not self.template:
            cv['resultsFolder'] += cv['experiment'] + '/'
        self.update_rep_result_folders(cv)

        # Standard labels:
        if cv['exp_type'] in self.STANDARD_LABELS:
            labels = self.STANDARD_LABELS[cv['exp_type']]
            cv['description'] = "The ENCODE "+labels['long']+" pipeline"
            cv['title'] = labels['long']
            cv['name'] = labels['short']
            if not self.no_refs:
                cv['name'] = labels['short'] + "_"+cv['genome']
                if not self.template:
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
                if not self.no_refs:
                    rep['title'] += " on "+cv['genome']
                    if not self.template:
                        rep['title'] += ", "+cv['gender']
                rep['name']  = cv['name']+"_"+rep['rep_tech']
            cv['title']   += " - "+cv['rep_tech']
            if not self.no_refs:
                cv['title']   += " on "+cv['genome']
                if not self.template:
                    cv['title']   += ", "+cv['gender']
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
        full_mapping = encd.get_full_mapping(exp_id,exp=exp)
        controls_expected = (self.CONTROL_FILE_GLOB != None)
        reps = encd.get_reps(exp_id, load_reads=True, exp=exp, full_mapping=full_mapping, \
                             control_locs=controls_expected)
        rep_techs = []
        if 'reps' in args and args.reps != None:
            for rep_tech in args.reps:
                # Normalize rep_tech string
                if '_' not in rep_tech:
                    rep_tech = rep_tech + '_1' # default to tech_rep 1
                if not rep_tech.startswith('rep'):
                    rep_tech = 'rep' + rep_tech
                rep_techs.append(rep_tech)
        elif args.br != None and args.br != 0:
            rep_tech = 'rep' + str(args.br) + '_' + str(args.tr)
            rep_techs.append(rep_tech)

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
                            print >> sys.stderr, "ERROR: Specify different replicates to compare (e.g. 'rep1_1 rep2_1')."
                            sys.exit(1)
                    print "ERROR: >> sys.stderr, requested '"+rep_tech+"' not found in ENCODE " + cv['experiment']
                    sys.exit(1)
                ltr = letters.popleft()
                cv_reps[ltr] = reps[found_ix]
                del reps[found_ix] # rep is only used once!
        else:
            # Assume combined if exactly 2 replicates are available
            for rep in reps:
                ltr = letters.popleft()
                cv_reps[ltr] = rep

        # Special case for srna ENCODE2:
        if cv['exp_type'] == "small-rna-seq":
            for ltr in cv_reps.keys():
                a_tailing = encd.has_a_tailing(self.exp,cv_reps[ltr]['rep_tech'])
                if a_tailing is not None:
                    cv_reps[ltr]['a_tailing'] = a_tailing

        # Assume combined reps if supported AND exactly 2 reps AND for different biological reps
        if self.PIPELINE_BRANCH_ORDER != None and 'COMBINED_REPS' in self.PIPELINE_BRANCH_ORDER:
            if len(cv_reps) == 2:
                if cv_reps['a']['br'] != cv_reps['b']['br']:
                    self.combined_reps = True
                elif self.compare_techreps and cv_reps['a']['br'] == cv_reps['b']['br'] \
                                           and cv_reps['a']['tr'] != cv_reps['b']['tr']:
                    self.combined_reps = True
        self.multi_rep = (len(reps) > 1)

        # A little more rep tidying
        for ltr in cv_reps.keys():
            rep = cv_reps[ltr]
            rep['branch_id'] = 'REP'  # MUST override for non-standard pipelines which have ['REP','COMBINED_REPS']
            rep['concat_id'] = 'reads'
            if 'paired_end' in rep and rep['paired_end']:
                rep['concat_id2'] = 'reads2'
            if self.detect_umi:
                (umi_found_true, barcodes) = encd.rep_is_umi(exp,rep)
                if umi_found_true:
                    rep['umi'] = 'yes'
                else:
                    rep['umi'] = 'no'
                if len(barcodes) > 0:
                    rep['barcode'] = barcodes[0]
                else:
                    rep['barcode'] = "undetected"

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


    def load_template_reps(self, args, cv):
        '''
        Creates template replicates building the cv['reps'] list from requested or all.
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

        # Determine rep_techs requested
        rep_techs = []
        if 'reps' in args and args.reps != None:
            for rep_tech in args.reps:
                # Normalize rep_tech string
                if '_' not in rep_tech:
                    rep_tech = rep_tech + '_1' # default to tech_rep 1
                if not rep_tech.startswith('rep'):
                    rep_tech = 'rep' + rep_tech
                rep_techs.append( (rep_tech,int(rep_tech[3:].split('_')[0]),int(rep_tech[3:].split('_')[1])) )
        elif args.br != None and args.br != 0:
            rep_tech = 'rep' + str(args.br) + '_' + str(args.tr)
            rep_techs.append( (rep_tech,int(args.br),args(args.tr)) )
        else:
            rep_techs = [ ('rep1_1',1,1), ('rep2_1',2,1) ]  # Default is 2 biological replicate pipeline run

        # Create rep dicts that will be used to build pipeline ins and outs
        cv_reps = {}
        letters = deque('abcdefghijklmnopqrstuvwxyz')
        processed = []
        for (rep_tech,br,tr) in rep_techs:
            if rep_tech in processed:
                print >> sys.stderr, "ERROR: requested %s twice." % rep_tech
                sys.exit(1)
            processed.append(rep_tech)
            rep = { 'br': br, 'tr': tr,'rep_tech': rep_tech }
            rep['organism'] = 'human'
            if args.genome == 'mm10':
                rep['organism'] = 'mouse'
            rep['sex'] = 'male' # Sorry, but the patriarchy has ruled that we must include chrY.
            rep['paired_end'] = args.pe  # Defaults to se
            rep['has_reads'] = False

            # Load read and control files, only if requested.
            run_type = None  # The files should be consistent as all single-end or paired-end, but some mapping got it wrong
            rep['fastqs'] = { "1": [], "2": [] }
            rep['controls'] = []
            rep['priors'] = {}
            rep['inputs'] = {}
            rep['inputs']['Reads1'] = []
            rep['inputs']['Reads2'] = []

            ltr = letters.popleft()
            cv_reps[ltr] = rep
            #print "* Added: cv_reps[%s] as %s" % (ltr, rep_tech)

        # Add combining replicates as appropriate
        # Assume combined reps if supported AND exactly 2 reps AND for different biological reps
        if self.PIPELINE_BRANCH_ORDER != None and 'COMBINED_REPS' in self.PIPELINE_BRANCH_ORDER:
            if len(cv_reps.keys()) == 2 and cv_reps['a']['br'] != cv_reps['b']['br']:
                self.combined_reps = True
        self.multi_rep = (len(cv_reps.keys()) > 1)

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
        for rep_id in sorted( cv['reps'].keys() ):  # Everything will be initiated
            rep = cv['reps'][rep_id]
            if "priors" not in rep:
                rep['priors'] = {}
                rep['inputs'] = {}


    def update_rep_result_folders(self, cv):
        ''' create input object for a step and extends the file_globs dict as appropriate.'''
        # NOT EXPECTED TO OVERRIDE
        for rep_id in sorted( cv['reps'].keys() ):
            rep = cv['reps'][rep_id]
            if self.template or rep_id == self.SEA_ID:  # SEA is the ultimate branch where all pipelines flow
                rep['resultsFolder'] = cv['resultsFolder']
            elif 'tributaries' not in 'rep':
                rep['resultsFolder'] = cv['resultsFolder'] + rep['rep_tech'] + '/'
            elif 'resultsFolder' in 'rep':
                rep['resultsFolder'] = cv['resultsFolder'] + rep['resultsFolder'] + '/'
            else:
                rep['resultsFolder'] = cv['resultsFolder']


    def build_applets_if_necessary(self):
        ''' create input object for a step and extends the file_globs dict as appropriate.'''
        # Build apps first (if necessary)...
        if self.build_apps:
            build_cmd = os.path.dirname(sys.argv[0]) + '/build_applets'
            # TODO: make sure to use the correct project self.proj_id
            if self.test:
                print "Requires building the applets with '%s' before this can succeed!" % build_cmd
            else:
                # TODO: Just request the applets being used.
                print "Requesting to build the applets with '%s' ... " % build_cmd
                subprocess.call([build_cmd])
                self.build_apps = False # Done with that!

    def build_a_step(self, applet, file_globs, proj_id):
        ''' create input object for a step and extends the file_globs dict as appropriate.'''
        # NOT EXPECTED TO OVERRIDE
        try:
            dxapp = dxpy.find_one_data_object(classname="applet", name=applet, project=proj_id,
                                              zero_ok=False, more_ok=False,describe=True)
        except IOError:
            print >> sys.stderr, "ERROR: Cannot find applet '"+applet+"' in "+proj_id
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
        self.build_applets_if_necessary()
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
            # NOTE: Useful for discovering pipeline json, which can be copied/modified for new launcher
            #pipe_steps, file_globs = self.build_simple_steps(pipe_path,self.proj_id,verbose=True)
            #sys.exit(1)

        for rep in self.psv['reps'].values():
            if rep["branch_id"] != branch_id:
                continue
            if type(pipe_path) == dict:
                assert 'se' in pipe_path and 'pe' in pipe_path
                if rep.get('paired_end',False):
                    rep['path'] = pipe_path['pe']
                else:
                    rep['path'] = pipe_path['se']
            else:
                rep['path'] = pipe_path
            # In a rare case, pruning of steps is done.
            if self.PRUNE_STEPS != None and len(self.PRUNE_STEPS) != 0:
                for step in self.PRUNE_STEPS:
                    if step in pipe_steps:
                        del pipe_steps[step]
                    if step in rep['path']:
                        rep['path'].remove(step)
            rep['steps'] = pipe_steps

        return [ pipe_steps, file_globs ]

    def find_prior_results(self,pipe_path,steps,results_folder,file_globs):
        '''Looks for all result files in the results folder.'''
        priors = {}
        for step in pipe_path:
            for file_token in steps[step]['results'].keys():
                fid = self.find_file(results_folder + file_globs[file_token],self.proj_id, \
                                                                                    recurse=False)
                if fid != None:
                    priors[file_token] = fid
        return priors


    def find_inputs_and_priors(self,branch_id,globs,verbose=False):
        '''Finds the inputs and priors for all reps for a given branch of the pipeline.'''
        # NOT EXPECTED TO OVERRIDE
        #verbose=True

        if not self.template:
            print "Checking for prior results..."
            # NOTE: priors is a dictionary of fileIds that will be used to determine stepsToDo
            #       and fill in inputs to workflow steps
            for rep in self.psv['reps'].values():
                if rep['branch_id'] != branch_id:
                    continue
                if not self.test:
                    if not dx.project_has_folder(self.project, rep['resultsFolder']):
                        self.project.new_folder(rep['resultsFolder'],parents=True)
                rep['priors'] = self.find_prior_results(rep['path'],rep['steps'],rep['resultsFolder'],globs)

        if not self.template:
            print "Checking for input files..."
            # Find all reads files and move into place
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
                rep['inputs']['Reads1'] = dx.find_and_copy_read_files(rep['priors'], \
                                                    rep['fastqs']['1'],  self.test, reads_token, \
                                                    rep['resultsFolder'], False, self.proj_id)
                # Note: rep['fastqs']['2'] and rep['inputs']['Reads2'] will be empty on single-end
                rep['inputs']['Reads2'] = dx.find_and_copy_read_files(rep['priors'], \
                                                    rep['fastqs']['2'], self.test, 'reads2', \
                                                    rep['resultsFolder'], False, self.proj_id)
                if verbose:
                    print "%s" % str(rep['fastqs']['1'])
                    print json.dumps(rep['inputs']['Reads1'],indent=4,sort_keys=True)
                    print "%s" % str(rep['fastqs']['2'])
                    print json.dumps(rep['inputs']['Reads2'],indent=4,sort_keys=True)

        for rep in self.psv['reps'].values():
            # If rep has tributaries, then look for those priors!
            if rep['branch_id'] != branch_id:
                continue
            if 'tributaries' in rep:
                self.find_combined_inputs(rep,rep['steps'],globs)


    def find_all_ref_files(self):
        '''Locates all reference files based upon organism and gender.'''
        # NOT EXPECTED TO OVERRIDE

        if not self.no_refs:
            print "Looking for reference files..."
            ref_priors = {}
            ref_priors = self.find_ref_files(ref_priors)

            # Add into all rep priors
            for rep in self.psv['reps'].values():
                rep['priors'].update(ref_priors)


    def find_all_app_non_file_outputs(self,verbose=False):
        '''Finds app non-file outputs, or the file parameters that match by name.'''
        # NOT EXPECTED TO OVERRIDE
        # Another obscure and difficult one!
        #verbose=True

        print "Checking for prior app non-file outputs..."

        # After inputs and priors have been discovered...
        # 1) walk through all reps
        for rep_id in sorted(self.psv["reps"].keys()):
            rep = self.psv["reps"][rep_id]
            if "params" not in rep:
                rep["params"] = {}
            # 2) walk through all steps in a rep
            for step_id in rep["steps"]:
                step = rep["steps"][step_id]
                # 3) if step has an "output_values" then for each output_value key:
                if "output_values" not in step:
                    continue
                if verbose:
                    print >> sys.stderr, "DEBUG: output_values in "+rep['rep_tech']+":"+step_id
                # Note: set of files and set of outputs, but any file could lead to job and job should have all outputs
                #       Still missing jobs or output not in job means look in files and only one file likely has an output
                job = None # only need to find once per step
                for file_token in step["results"]:
                    # 4) if a result file is found...
                    if file_token not in rep["priors"]:
                        continue
                    fid = rep["priors"][file_token]
                    for out_key in step["output_values"].keys():
                        if out_key in rep["params"]:
                            continue
                        out_name = step["output_values"][out_key]
                        out_value = None
                        # 5) look up the job, then look up the output value matching this one
                        if job == None:
                            job = dx.job_from_fid(fid)
                        if job:
                            if "output" in job and out_name in job["output"]:
                                out_value = job["output"][out_name]
                                if verbose:
                                    print >> sys.stderr, "DEBUG: found in job '"+out_name+"' = "+str(out_value)+", saved as "+ \
                                                                                           rep_id+"['params']['"+out_key+"']"
                        # 6) if not found then check in this file
                        if out_value == None:
                            out_value = dx.file_get_property(out_name,fid)
                            if out_value != None:
                                if verbose:
                                    print >> sys.stderr, "DEBUG: found in file '"+out_name+"' = "+str(out_value)+", saved as "+ \
                                                                                           rep_id+"['params']['"+out_key+"']"
                        # 7) If found add value to rep["params"][out_key] and expect there are no name conflicts
                        if out_value != None:  # Put this thing somewhere!
                            rep["params"][out_key] = out_value


    def determine_steps_to_run(self,pipe_path, steps, priors, deprecate, force=False, verbose=False):
        '''Determine what steps need to be done, base upon prior results.'''
        # NOT EXPECTED TO OVERRIDE
        #verbose=True  # Very useful when verifying new/updated launcher
        will_create = []
        steps_to_run = []
        for step in pipe_path:
            if self.template:
                steps_to_run += [ step ]
                if verbose:
                    print "- Adding step '"+step+"' because all steps are needed when templating."
                continue
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
                    continue
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
                    if inp in priors and priors[inp] in self.deprecate: # This should catch deprecations across rep boundaries.
                        steps_to_run += [ step ]
                        if verbose:
                            print "- Adding step '"+step+"' due to tributary results being deprecated."
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
                        self.deprecate.append(priors[result]) # needed to ensure tributaries affect river reps.
                        if verbose:
                            print "  - Will deprecate '"+dx.file_path_from_fid(priors[result])+"'\."
                        del priors[result]
                        # if results are in folder, then duplicate files cause a problem!
                        # So add to 'deprecate' to move or remove before launching

        # Now make sure the steps can be found, and error out if not.
        self.build_applets_if_necessary()
        for step in steps_to_run:
            app = steps[step].get('app',step)
            dxApp = dxpy.find_data_objects(classname='file', name=app, name_mode='exact',
                                                        project=self.proj_id, return_handler=False)
            if dxApp == None:
                print >> sys.stderr, "ERROR: failure to locate app '"+app+"'!"
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
        '''If combined_reps step, look for input from results of replicate level step.'''
        # NOT EXPECTED TO OVERRIDE

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
            # Expectation that when combining results, dx app's input names either end in '_A', '_ABC', etc. or '_set'
            results_array = []
        elif inp_token[-2] == '_' and inp_token.lower()[-1] in letters:
            looking_for = inp_token.lower()[-1]
        elif looking_for == None:
            # NOT a warning: no current means to distinguish between expecting a prior from not expecting one
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


    def find_mystery_param(self,rep,rep_id,step_id,param_name,link_later=False,verbose=False):
        '''Finds parameters that are outputs of previous steps.'''
        # NOT EXPECTED TO OVERRIDE, though only needed for dnaseLaunch.py so far

        # Mystery params come from the non-file outputs of other steps (even from other reps!).
        # If a step definition contains an "output_values" parameter then:
        # - find_all_app_non_file_outputs will look for already run step values and place them in rep['params']
        # - create_or_extend_workflow will create links in rep['prevStepResults'] for steps not yet run
        # If this step contains "param_links", then each param_spec should say what rep to look for: 'self', 'sister'
        #verbose=True

        # Not easy, that is for sure! "param_links": { "target_size": {"name":"reads_filtered", "rep":"sister"} }
        if "param_links" not in rep["steps"][step_id]:
            return None
        param_links = rep["steps"][step_id]["param_links"]
        if param_name not in param_links:
            return None
        param_spec = param_links[param_name]
        mystery_rep_id = None
        if param_spec["rep"] == "sister": # Hellena's 'brother sistra'?  Catch the reference anyone?
            # How to find "brother sistra"?
            if "sister" not in rep:
                print "ERROR: Looking for mystery parameter but can't find sister of rep '"+rep_id+"'."
                sys.exit(0)
            mystery_rep_id = rep["sister"]
        elif param_spec["rep"] == "self":
            mystery_rep_id = rep_id
        else: # or "parent" or "aunt" ??
            print "ERROR: Looking for mystery parameter but don't (yet) support rep '"+param_spec["rep"]+"'."
            sys.exit(0)

        mystery_rep = self.psv["reps"][mystery_rep_id]
        if verbose:
            print >> sys.stderr, "DEBUG: Found rep "+rep['rep_tech']+"'s "+param_spec["rep"]+" as mystery rep '"+mystery_rep['rep_tech']+"'"
        mystery_value = None
        if "prevStepResults" in mystery_rep and param_spec["name"] in mystery_rep["prevStepResults"]:
            mystery_value = mystery_rep["prevStepResults"][param_spec["name"]]  # not yet run result (dx link)
            if verbose:
                print >> sys.stderr, "DEBUG: Found mystery param '"+param_name+"' as '"+param_spec["name"]+\
                      "' from mystery rep '"+mystery_rep['rep_tech']+"' previous results, with value: " + str(mystery_value)
            return mystery_value
        elif "params" in mystery_rep and param_spec["name"] in mystery_rep["params"]:
            mystery_value = mystery_rep["params"][param_spec["name"]] # prior run result (actual string value)
            if verbose:
                print >> sys.stderr, "DEBUG: Found mystery param '"+param_name+"' as '"+param_spec["name"]+\
                      "' from mystery rep '"+mystery_rep['rep_tech']+"' with value: " + str(mystery_value)
            return mystery_value

        if not link_later or param_spec["rep"] == "self":  # Linking to self must be to prior step in rep
            return None

        # Avoid failure to link BOTH mystery 'sister' params because one is always prior to the other
        # Solution is to put in a known token '@link_later@', then have a final pass after building whole workflow
        # This link_later strategy should only be used if the mystery_step that produces the value can be found.
        mystery_out = None
        for mystery_step_key in mystery_rep["steps"].keys():
            mystery_step = mystery_rep["steps"][mystery_step_key]
            if "output_values" in mystery_step and param_spec["name"] in mystery_step["output_values"]:
                mystery_out = mystery_step["output_values"][param_spec["name"]]
                if verbose:
                    print >> sys.stderr, "DEBUG: Found mystery rep '"+mystery_rep_id+"' and step "+mystery_step_key+"' has '"+ \
                                                                        param_name+" as '"+mystery_out+"'"
                break
        if mystery_out != None:
            mystery_value = '@link_later@'
            if verbose:
                print >> sys.stderr, "DEBUG: Found mystery param '"+param_name+"' as '"+param_spec["name"]+\
                      "' from rep '"+param_spec["rep"]+"' with value: " + str(mystery_value)
            return mystery_value

        return None


    def wf_find_file_input(self,rep,step_id,file_token,expect_set,verbose=False):
        '''When buidling a workflow, finds a dx defined step file input (which may be for a set)
           from prior step results, tributary results, or prior files found.'''
        # NOT EXPECTED TO OVERRIDE
        #verbose=True

        # Now try to find prior results or existing files to fill in this input
        file_input = []
        expect_count = 1
        if 'tributaries' in rep:
            # First: check in previous branch (e.g. rep level pipeline) is input to the combining step
            if expect_set and \
               (file_token.endswith('_A') or file_token.endswith('_B') or file_token.endswith('_ABC')):
                expect_count = len(rep['tributaries'])
            prev_input = self.find_results_from_prev_branch(rep, file_token, expect_set)
            if prev_input != None:
                if isinstance(prev_input, list):
                    file_input.extend( prev_input )
                else:
                    file_input.append( prev_input )

            if verbose:
                print >> sys.stderr, "DEBUG: Found %d file_inputs from prev_branch" % len(file_input)

        if len(file_input) != expect_count and file_token in rep['prevStepResults']:
            if isinstance(rep['prevStepResults'][file_token], list):
                file_input.extend( rep['prevStepResults'][file_token] )
            else:
                file_input.append( rep['prevStepResults'][file_token] )
            if verbose:
                print >> sys.stderr, "DEBUG: Now %d file_inputs after prevStepResults" % len(file_input)

        if len(file_input) != expect_count:
            alt_token = file_token
            if not self.template:
                if alt_token not in rep['priors'] and expect_set:
                    alt_token = file_token + "_set"  # Sometimes sets expect  the token to end in '_set'!
                if alt_token not in rep['priors'] and file_token.startswith('reads'):
                    alt_token = "reads"  # Sometimes sets are just reads.
                if alt_token not in rep['priors'] and alt_token == "reads" and expect_set:
                    alt_token = "reads_set"  # Sometimes reads are reads_sets. # FIXME: should really drop "_set" logic
            if alt_token == "reads" and file_token == "reads2" and not rep['paired_end']:
                return None # This special case can come when applet is built for PE and SE reads with optional reads2 param
            #print >> sys.stderr, "DEBUG: alt_token '%s' file_token '%s'" % (alt_token,file_token)
            if alt_token in rep['priors']:
                if isinstance(rep['priors'][alt_token], list):
                     if expect_set:
                        for fid in rep['priors'][alt_token]:
                            file_input.append( dx.FILES[fid] )
                     ##else:
                     #   print >> sys.stderr, "ERROR: Not expecting input set and found "+str(len(rep['priors'][alt_token]))+" file(s)."
                     #   print >> sys.stderr, json.dumps(rep['priors'],indent=4,sort_keys=True)
                     #   sys.exit(1)
                else:
                    file_input.append( dx.FILES[ rep['priors'][alt_token] ] )
                if verbose:
                    print >> sys.stderr, "DEBUG: Now %d file_inputs after priors, with '%s'" % (len(file_input),file_token)

        if not self.template:
            if len(file_input) == 0:
                print >> sys.stderr, "ERROR: step '%s'  can't find input '%s' paired:%s" % (step_id,file_token,str(rep['paired_end']))
                #print >> sys.stderr, json.dumps(rep['priors'],indent=4,sort_keys=True)
                sys.exit(1)
            if 'tributaries' in rep and len(file_input) != expect_count:
                print >> sys.stderr, "ERROR: step '"+step_id+"' input '"+file_token+"' expects %d files but found %d." % \
                                                                                    (expect_count,len(file_input))
                #print >> sys.stderr, json.dumps(rep['priors'],indent=4,sort_keys=True)
                sys.exit(1)

        if len(file_input) == 0:
            return None
        if not expect_set:
            return file_input[0]
        return file_input


    def wf_find_param_input(self,rep,rep_id,step_id,param,app_param,verbose=False):
        '''When buidling a workflow, finds a dx defined step param input (sets not yet supported)
           from rep, psv or in rare cases another step'''
        # NOT EXPECTED TO OVERRIDE

        # Now try to find prior results or existing files to fill in this input
        param_input = None
        if param in rep:
            param_input = rep[param]
        #elif "params" in rep and param in rep["params"]: # watch it... this may have cross-fertilizing params
        #    param_input = rep["params"][param]
        elif param in self.psv:
            param_input = self.psv[param]
        else:
            # Rare case of non-file outputs of one step being inputs of another
            param_input = self.find_mystery_param(rep, rep_id, step_id, app_param,link_later=True)
        if param_input == None and not self.template:
            print >> sys.stderr, "ERROR: step '"+step_id+"' unable to locate '"+param+ \
                                              "' in pipeline specific variables (psv)."
            #print >> sys.stderr, json.dumps(rep,indent=4)
            sys.exit(1)

        return param_input


    def create_or_extend_workflow(self,rep, rep_id, wf=None,app_proj_id=None):
        '''
        This function will populate or extend a workflow with the steps in rep['stepsToDo'] and return
        the workflow unlaunched.   It relies on the rep object containing *almost* all information to run
        the PIPELINE_BRANCH associated with the rep.  The rep should contain rep['priors'] defining existing files;
        rep['steps'] objects defining applet, input and output requirements; etc.  A few common parameters will
        be pulled from the psv (pipeline specific variables) object.
        '''
        # NOT EXPECTED TO OVERRIDE

        # Make sure results folder exists first.
        if not self.test:
            if not dx.project_has_folder(self.project, rep['resultsFolder']):
                self.project.new_folder(rep['resultsFolder'],parents=True)
        self.build_applets_if_necessary()

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
        if not self.test and wf == None:
            wf = dxpy.new_dxworkflow(title=wf_name,name=wf_name,folder=wf_folder,
                                           project=self.proj_id,description=self.psv['description'])

        # NOTE: prevStepResults dict contains links to result files to be generated by previous steps
        if 'prevStepResults' not in rep:
            rep['prevStepResults'] = {}
        steps = rep['steps']
        for step in rep['stepsToDo']:
            step_link_later = False
            app_name = steps[step].get('app',step)
            app = dx.find_applet_by_name(app_name, app_proj_id)
            inp_defs = app.describe().get('inputSpec') or []
            app_inputs = {}

            # file inputs
            for file_token in steps[step]['inputs'].keys():
                app_inp = steps[step]['inputs'][file_token]
                if app_inp == "reads2" and not rep['paired_end']:
                    continue
                # Need inp_def to see if set is expected
                for inp_def in inp_defs:
                    if inp_def["name"] == app_inp:
                        break
                if inp_def["name"] != app_inp:
                    print >> sys.stderr, "ERROR: Pipeline definition for applet '"+app_name+"' file_token '"+file_token+"'"
                    print >> sys.stderr, "value of %s not found in DX definition." % app_inp
                    sys.exit(1)
                expect_set = (inp_def["class"] == 'array:file')
                inp_file = self.wf_find_file_input(rep,step,file_token,expect_set)
                if inp_file != None:
                    app_inputs[ app_inp ] = inp_file

            # Non-file app param inputs
            if 'params' in steps[step]:
                for param in steps[step]['params'].keys():
                    app_param = steps[step]['params'][param]
                    inp_param = self.wf_find_param_input(rep,rep_id,step,param,app_param)
                    if inp_param != None:
                        app_inputs[ app_param ] = inp_param
                        # Since a mystery parameter could be from a sister rep, it may require linking later
                        if isinstance(app_inputs[ app_param ],str) and app_inputs[ app_param ] == '@link_later@':
                            step_link_later = True  # Need to save this step's stage for later
                            # Need a final pass of workflow once all branches are built
                            if self.link_later == None:
                                self.link_later = []
                            later_link = { 'rep': rep, 'rep_id': rep_id, 'step_id': step, 'param': app_param }
                            self.link_later.append(later_link)

            # Now we are ready to add wf stage
            if not self.test:
                if self.template:
                    stage_id = wf.add_stage(app, stage_input=app_inputs) # Templates should not fill in a result folder!
                else:
                    stage_id = wf.add_stage(app, stage_input=app_inputs, folder=rep['resultsFolder'])
                if step_link_later:  # At least one parameter will requre linking later so save stage_id
                    steps[step]['stage_id'] = stage_id

            # outputs, which we will need to link to
            for file_token in steps[step]['results'].keys():
                app_out = steps[step]['results'][file_token]
                if self.test:
                    rep['prevStepResults'][ file_token ] = 'fake-for-testing'
                else:
                    rep['prevStepResults'][ file_token ] = dxpy.dxlink({'stage': stage_id, \
                                                                'outputField': app_out })
            # need to add output_values to previous results also (to support param_link requests).
            if "output_values" in steps[step]:
                for out_token in steps[step]["output_values"].keys():
                    app_out = steps[step]["output_values"][out_token]
                    if self.test:
                        rep['prevStepResults'][ out_token ] = 'fake-for-testing'
                    else:
                        rep['prevStepResults'][ out_token ] = dxpy.dxlink({'stage': stage_id, \
                                                                     'outputField': app_out })

        return wf


    def workflow_final_pass(self, wf, app_proj_id=None,verbose=False):
        '''
        Final pass of workflow just in case any '@link_later@' parameters need to be filled in.
        '''
        # NOT EXPECTED TO OVERRIDE
        if self.link_later == None:
            return

        # Since non-file outputs from one step might be input params to another step and since
        # those outputs might be in 'sister' branches of the workflow, this final pass may be
        # needed to fill in the links that were not available when the wf was built branch by branch.
        # TODO: This logic could also be used for file results->inputs, but is not currently.
        #verbose=True

        print "Final pass through workflow..."

        for later_link in self.link_later:
            rep       = later_link['rep']
            rep_id    = later_link['rep_id']
            step_id   = later_link['step_id']
            app_param = later_link['param']
            if verbose:
                print >> sys.stderr, "DEBUG: looking for '"+rep['rep_tech']+"' step '"+step_id+"' param '"+app_param+"'"
            mystery_param = self.find_mystery_param(rep,rep_id,step_id, app_param)
            if mystery_param == None:
                print >> sys.stderr, "ERROR: Link later for '"+rep['rep_tech']+"' step '"+step_id+"' param '"+app_param+"' failed to find link."
                sys.exit(1)
            # Could also assert that mystery_param is a dx link
            steps = rep['steps']
            step = steps[step_id]
            stage_id = step['stage_id'] # should not be a key error!
            stage = wf.get_stage(stage_id)
            if verbose:
                print "Stage:"
                print stage
            app_inputs = stage['input']
            if verbose:
                print >> sys.stderr, "DEBUG: found '"+app_inputs[app_param]+"' which is being replaced with '"+mystery_param+"'"
            app_inputs[app_param] = mystery_param
            wf.update_stage(stage_id, stage_input=app_inputs)

        return


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
        if self.multi_rep or self.combine_one_or_more:
            for ltr in sorted( self.psv['reps'].keys() ):
                rep = self.psv['reps'][ltr]
                rep_tech_msg = rep['rep_tech']
                for input_type in sorted( rep['inputs'].keys() ):
                    input_type_msg = input_type
                    if len(rep['inputs'][input_type]) > 0:
                        for fid in rep['inputs'][input_type]:
                            print "  - "+rep_tech_msg+' '+input_type_msg+' '+dx.file_path_from_fid(fid)
                            rep_tech_msg = '      '  # Trick to make output less cluttered
                            input_type_msg = '      '
        else:
            rep_tech_msg = run['rep_tech']
            for input_type in sorted( run['inputs'].keys() ):
                input_type_msg = input_type
                if len(run['inputs'][input_type]) > 0:
                    for fid in run['inputs'][input_type]:
                        print "  - "+rep_tech_msg+' '+input_type_msg+' '+dx.file_path_from_fid(fid)
                        rep_tech_msg = '      '
                        input_type_msg = '      '

        if not self.no_refs:
            print "- Reference files:"
            # Should be enough to show the ref priors from 'a' rep for single-rep, multi or combined
            for token in self.psv['ref_files']:
                print "  " + dx.file_path_from_fid(self.psv['reps']['a']['priors'][token],True)

        print "- Results written to: " + self.psv['project'] + ":" +run['resultsFolder']
        if not self.template:
            if self.multi_rep:
                for rep_id in sorted( self.psv['reps'].keys() ):
                    if rep_id != self.SEA_ID:
                        rep = self.psv['reps'][rep_id]
                        if rep['branch_id'] not in ["REP","TECH_REP"]:
                            print "                      " + self.psv['project'] + ":" +rep['resultsFolder'] + " (combining)"
                        else:
                            print "                      " + self.psv['project'] + ":" +rep['resultsFolder']

        print "- Steps to run:"
        to_run_count = 0
        if self.multi_rep or self.combine_one_or_more:
            for ltr in sorted( self.psv['reps'].keys() ):
                rep = self.psv['reps'][ltr]
                for step in rep['path']:
                    if step in rep['stepsToDo']:
                        print "  * "+rep['rep_tech'] +': '+ rep['steps'][step].get('app',step)+" will be run"
                    else:
                        if not step.find('concat') == 0:
                            print "    "+rep['rep_tech'] +': '+ rep['steps'][step].get('app',step)+" has already been run"
                to_run_count += len(rep['stepsToDo'])
        else:
            to_run_count += len(run['stepsToDo'])
            for step in run['path']:
                if step in run['stepsToDo']:
                    print "  * "+run['steps'][step].get('app',step)+" will be run"
                else:
                    if not step.find('concat') == 0:
                        print "    "+run['steps'][step].get('app',step)+" has already been run"
        if to_run_count == 0:
            print "* All expected results are in the results folder, so there is nothing to do."
            print "  If this experiment/replicate needs to be rerun, then use the --force flag to "
            print "  rerun all steps; or remove suspect results from the folder before relaunching."
            sys.exit(0)

        if not self.template:
            if self.multi_rep or self.combine_one_or_more:
                for ltr in sorted( self.psv['reps'].keys() ):
                    rep = self.psv['reps'][ltr]
                    if len(rep['deprecate']) > 0:
                        deprecated = rep['resultsFolder']+"/deprecated/"
                        if rep['resultsFolder'].endswith('/'):
                            deprecated = run['resultsFolder']+"deprecated/"
                        print "Will move "+str(len(rep['deprecate']))+" prior result file(s) to '" + deprecated+"'."
                        for fid in rep['deprecate']:
                            print "  " + dx.file_path_from_fid(fid)
            else:
                if len(run['deprecate']) > 0:
                    deprecated = run['resultsFolder']+"/deprecated/"
                    if run['resultsFolder'].endswith('/'):
                        deprecated = run['resultsFolder']+"deprecated/"
                    print "Will move "+str(len(run['deprecate']))+" prior result file(s) to '" + deprecated+"'."
                    for fid in run['deprecate']:
                        print "  " + dx.file_path_from_fid(fid)


    def workflow_report_and_build(self,run,proj_id=None,verbose=False):
        '''Builds the multi-rep/combined-rep workflow from parts, reporting plans and returning wf.'''
        # NOT EXPECTED TO OVERRIDE
        # Report the plans
        if self.multi_rep:
            self.report_run_plans()
        else:
            self.report_run_plans(run)

        if not self.template:
            print "Checking for currently running analyses..."
            self.check_run_log(run['resultsFolder'], proj_id, verbose=True)

            # Move old files out of the way...
            if self.multi_rep or self.combine_one_or_more:
                for rep_id in sorted( self.psv['reps'].keys() ):
                    rep = self.psv['reps'][rep_id]
                    if len(rep['deprecate']) > 0 and not self.test:
                        deprecated = rep['resultsFolder']+"deprecated/"
                        print "Moving "+str(len(rep['deprecate']))+" "+rep['rep_tech']+" prior result file(s) to '"+ \
                                                                                            deprecated+"'..."
                        dx.move_files(rep['deprecate'],deprecated,proj_id)
            if not self.multi_rep and not self.combine_one_or_more:
                if len(run['deprecate']) > 0 and not self.test:
                    deprecated = run['resultsFolder']+"deprecated/"
                    print "Moving "+str(len(run['deprecate']))+" "+run['rep_tech']+" prior result file(s) to '"+ \
                                                                                        deprecated+"'..."
                    dx.move_files(run['deprecate'],deprecated,proj_id)

        # Build the workflow...
        wf = None
        if self.multi_rep or self.combine_one_or_more:
            for rep_id in sorted( self.psv['reps'].keys() ):
                rep = self.psv['reps'][rep_id]
                dotdotdot = "..."
                if rep['branch_id'] not in ["REP","TECH_REP"]:
                    dotdotdot = " (combining)..."
                if len(rep['stepsToDo']) > 0: # Going to get noisy with multiple reps
                    if self.test:
                        print "Testing workflow assembly for "+rep['rep_tech']+dotdotdot
                    else:
                        print "Assembling workflow for "+rep['rep_tech']+dotdotdot
                    wf = self.create_or_extend_workflow(rep, rep_id, wf=wf)
        if not self.multi_rep and not self.combine_one_or_more:
            if len(run['stepsToDo']) > 0:
                if self.test:
                    print "Testing workflow assembly for "+run['rep_tech']+"..."
                else:
                    print "Assembling workflow for "+run['rep_tech']+"..."
                wf = self.create_or_extend_workflow(run, None, wf=wf)

        # Need a final pass over the whole workflow to patch in params set to '@link_later@'
        if wf != None:  # Only non-test case
            self.workflow_final_pass(wf)

        return wf


    def check_run_log(self,results_folder,proj_id,verbose=False):
        '''Checks for currently running jobs and will exit if found.'''
        # NOT EXPECTED TO OVERRIDE
        run_log_path = results_folder + '/' + dx.RUNS_LAUNCHED_FILE
        log_fids = self.find_file(run_log_path,proj_id,multiple=True,recurse=False)
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
                            msg="ERROR: Exiting: Can't launch because prior run ["+run_id[0]+"] "
                            if len(run_id) > 1:
                                msg+="("+run_id[1]+") "
                            msg+= "has not finished (currently '"+state+"')."
                            print >> sys.stderr, msg
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
        run_log_path = results_folder + '/' + dx.RUNS_LAUNCHED_FILE
        old_fids = self.find_file(run_log_path,self.proj_id,multiple=True,recurse=False)
        new_fh = dxpy.new_dxfile('w',project=self.proj_id,folder=results_folder, \
                                                                  name=dx.RUNS_LAUNCHED_FILE)
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
            print >> sys.stderr, "ERROR: failure to assemble workflow!"
            sys.exit(1)

        if self.template:
            print "Template workflow '" + wf.name + "' has been assembled in "+run['resultsFolder'] + "."
        elif ignition:
            print "Launch sequence initiating..."
            wf_run = wf.run({}, project=self.proj_id,priority="normal")
            if wf_run == None:
                print >> sys.stderr, "ERROR: failure to lift off!"
                sys.exit(1)
            else:
                print "  We have liftoff!"
                wf_dict = wf_run.describe()
                self.log_this_run(wf_dict['id'],run['resultsFolder'])
                print "  Launched " + wf_dict['id']+" as '"+wf.name+"'"
                encd.exp_patch_internal_status(self.psv['experiment'], 'processing')

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
            self.find_inputs_and_priors(branch_id,file_globs)

        # finding experiment specific control files in a stadardized way
        self.find_all_control_files()

        # finding pipeline specific reference files in a stadardized way
        self.find_all_ref_files()

        # Look for any prior app outputs that may be used in later steps
        self.find_all_app_non_file_outputs()

        # deternine steps to run in a stadardized way
        self.determine_steps_needed(args.force)

        # Preperation is done. Now build up multi-rep workflow
        if self.multi_rep or self.combined_reps:
            run = self.psv
        else:
            run = self.psv['reps']['a']
        wf = self.workflow_report_and_build(run,self.proj_id)

        # Exit if test only
        if self.test:
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

