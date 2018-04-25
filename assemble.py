#!/usr/bin/env python2.7
# assemble.py 0.0.2

import argparse,os, sys, json
import subprocess
#import urlparse, itertools, logging
#from datetime import datetime
#import Lib/fnmatch.py
import fnmatch

import dxpy
#import dxencode
import dx
import encd

### TODO:
#   1) Assemble *could* create a workflow_run object that represents a request for launch.
#   2) Ideally, a workflow_run object could be created on encoded as a request for assemble/launch.

# Assemble is meant to be run directly for most experiment types.  That is, no derived class should be needed.
#   Its purpose is to recognize all files (for a given experiment) that are available in ENCODEd and are needed to
#   complete a pipeline run (including fastqs and intermediate results); then if they are not already on DX copy them there.
#   Ideally assemble.py is run for a single experiment, but can act on a list of experiments.  The files are "fetched" to
#   DX by a DX applet and assemble.py waits for the fetch to complete for one experiment, before moving on to the next.
#   After a successful fetch, assemble.py can kick-off (ignite) the appropriate Launcher, so that the pipeline will be run.

class Assemble(object):
    '''
    Assembles all input files in dx for a given accessioned experiment and replicate (or combined replicates).
    '''

    SERVER_DEFAULT = 'www'
    '''At this time there is no need to use the any but the one true server for assembling.'''

    FOLDER_DEFAULT = '/runs/'
    '''This the default location for creating experiment folders on dnanexus.'''

    FILE_STATUSES_ACCEPTED = [ "released", "in progress" ]
    '''By default only 'released' files may be assembled.  Use --status-accepted to allow others.'''

    EXPERIMENT_TYPES_SUPPORTED = [ 'long-rna-seq', 'small-rna-seq', 'rampage', 'dnase-seq', "dna-me" ] #,"chip-seq" ]
    '''This module supports only these experiment (pipeline) types.'''

    LAUNCHERS = {   'long-rna-seq':     '../long-rna-seq-pipeline/dnanexus/lrnaLaunch.py',
                    'small-rna-seq':    '../long-rna-seq-pipeline/dnanexus/small-rna/srnaLaunch.py',
                    'rampage':          '../long-rna-seq-pipeline/dnanexus/rampage/rampageLaunch.py',
                    'dnase-seq':        '../dnase_pipeline/dnanexus/dnaseLaunch.py',
                    'dna-me':           '../dna-me-pipeline/dnanexus/dmeLaunch.py',
                }

    # Pipeline files includes inputs and results.  To assemble the files, there is no need to understand step order
    # dependencies.  Even replicate vs. combined is not important information as folder destinations can be discovered
    # from submitted_file_name paths.
    # For each "output_type" there will be one or more globs to recognize files.  Input files really don't need globs.
    PIPELINE_FILES = {
         "long-rna-seq": {
            "inputs":  { "reads": [ "*.fastqs.gz" ] },
            "results": { "alignments":                          [ "*_star_genome.bam",         "*_tophat.bam"          ],
                         "signal of all reads":                 [ "*_star_genome_all.bw",      "*_tophat_all.bw"       ],
                         "signal of unique reads":              [ "*_star_genome_uniq.bw",     "*_tophat_uniq.bw"      ],
                         "minus strand signal of all reads":    [ "*_star_genome_minusAll.bw", "*_tophat_minusAll.bw"  ],
                         "plus strand signal of all reads":     [ "*_star_genome_plusAll.bw",  "*_tophat_plusAll.bw"   ],
                         "minus strand signal of unique reads": [ "*_star_genome_minusUniq.bw","*_tophat_minusUniq.bw" ],
                         "plus strand signal of unique reads":  [ "*_star_genome_plusUniq.bw", "*_tophat_plusUniq.bw"  ],
                         "transcriptome alignments":   [ "*_star_anno.bam"          ],
                         "gene quantifications":       [ "*_rsem.genes.results"     ],
                         "transcript quantifications": [ "*_rsem.isoforms.results"  ] }
        },
        "small-rna-seq": {
            "inputs":  { "reads": [ "*.fastqs.gz" ] },
            "results": { "alignments":                 [ "*_star_genome.bam"    ],
                         "multi-read minus signal":    [ "*_small_minusAll.bw"  ],
                         "multi-read plus signal":     [ "*_small_plusAll.bw"   ],
                         "unique minus signal":        [ "*_small_minusUniq.bw" ],
                         "unique plus signal":         [ "*_small_plusUniq.bw"  ] }
        },
        "rampage": {
            "inputs":  { "reads": [ "*.fastqs.gz" ] },
            "results": { "alignments":                 [ "*_star_marked.bam" ],
                         "minus strand signal of all reads":    [ "*_5p_minusAll.bw"  ],
                         "plus strand signal of all reads":     [ "*_5p_plusAll.bw"   ],
                         "minus strand signal of unique reads": [ "*_5p_minusUniq.bw" ],
                         "plus strand signal of unique reads":  [ "*_5p_plusUniq.bw"  ],
                         "signal of all reads":                 [ "*_5p_all.bw" ],
                         "signal of unique reads":              [ "*_5p_uniq.bw"  ],
                         "transcription start sites":           [ "*_peaks.bed.gz", "*_peaks.gff.gz", "*_peaks.bb"  ],
                         "gene quantifications":                [ "*_quant.tsv"          ] }
            # TODO: "controls"?  Punt at this time
            # TODO: "references" Probably not.  Assume that they are in place.
        },
        "dnase-seq": { # TODO: Flesh out the DNase results
            "inputs":  { "reads": [ "*.fastqs.gz" ] },
            "results": { "unfiltered alignments":  [ "*_bwa_techrep.bam"    ] }
        },
        "dna-me": { # TODO: Flesh out the DNase results
            "inputs":  { "reads": [ "*.fastqs.gz" ] },
            "results": { "alignments":                 [ "*_bismark_techrep.bam" ] }
        },
    }

    def __init__(self):
        '''
        Assemble expects to be the only class for assembling experiments on for pipeline types.
        '''
        self.args = {} # run time arguments
        self.server_key = self.SERVER_DEFAULT  # TODO: replace with self.encd.server_key when Encd class is created
        self.proj_name = None
        self.project = None
        self.proj_id = None
        self.genome = None
        self.annotation = None
        self.exp = {}  # Will hold the encoded exp json
        self.full_mapping = None
        self.psv = {} # will hold pipeline specific variables.
        self.statuses_accepted = self.FILE_STATUSES_ACCEPTED
        print # TEMPORARY: adds a newline to "while retrieving session configuration" unknown error

    def get_args(self,parse=True):
        '''Parse the input arguments.'''
        # Change this description if there is a derived version of Assemble
        ap = argparse.ArgumentParser(description="Handles assembly of files needed to launched pipeline runs " +
                    "for supported experiment types. Can be run repeatedly and will only try to " +
                    "copy the files from ENCODEd that are not already on dnanexus in the expected loaction, All files " +
                    "are assembled to folder /<resultsLoc>/<experiment> and replicate sub-folders named as " +
                    "<experiment>/rep<biological-replicate>_<technical-replicate>.")

        ap.add_argument('-e', '--experiments',
                        help='One or more ENCODED experiment accessions',
                        nargs='+',
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

        ap.add_argument('-i','--inputs_only',
                        help='Only copy input files (e.g. fastqs) to dnanexus',
                        action='store_true',
                        required=False)

        # Easier to make it whole experiment or experiment/replicate
        #ap.add_argument('-cr','--cr','--combine-replicates',
        #                help="Combine or compare two replicates (e.g.'1 2_2').'",
        #                nargs='+',
        #                required=False)

        ap.add_argument('--project',
                        help="Project to run analysis in (default: '" + \
                                                        dx.env_get_current_project() + "')",
                        required=False)

        # Don't anticpate dealing with reference file
        #ap.add_argument('--refLoc',
        #                help="The location to find reference files (default: '" + \
        #                    dx.REF_PROJECT_DEFAULT + ":" + dx.REF_FOLDER_DEFAULT + "')",
        #                default=dx.REF_FOLDER_DEFAULT,
        #                required=False)

        ap.add_argument('-g','--genome',
                        help="The genome assembly for folder organization (default: GRCh38)",
                        default=None,
                        required=False)

        ap.add_argument('-a','--annotation',
                        help='Optionally enter annotation to help organize folders',
                        default=None,
                        required=False)

        ap.add_argument('-f','--folder',
                        help="The location to place experiment folders (default: " + \
                                                "'<project>:" + self.FOLDER_DEFAULT + "')",
                        default=self.FOLDER_DEFAULT,
                        required=False)

        # Like splashdown either test or run.  No third option!
        #ap.add_argument('--run',
        #                help='Run the workflow after assembling it.',
        #                action='store_true',
        #                required=False)

        ap.add_argument('--server',
                        help="Server to download files from (default: '" + self.SERVER_DEFAULT + "')",
                        default=self.SERVER_DEFAULT,
                        required=False)

        ap.add_argument('--sa','--status-accepted',
                        help='Optionally allow additional file statuses.',
                        nargs='+',
                        default=None,
                        required=False)

        ap.add_argument('-l','--launch',
                        help='Ignite the launcher after files have been assembled.',
                        action='store_true',
                        required=False)

        ap.add_argument('--test',
                        help='Test run only, do not assemble anything.',
                        action='store_true',
                        required=False)

        if parse:
            return ap.parse_args()
        else:
            return ap

    def load_variables(self,args,key='default'):
        '''Loads common variables to self.'''
        self.test = args.test
        self.inputs_only = args.inputs_only
        self.server_key = args.server
        encd.set_server_key(self.server_key) # TODO: change to self.encd = Encd(self.server_key)
        if self.server_key != "test":
            self.acc_prefix = "ENCFF"
        self.proj_name = dx.env_get_current_project()
        if self.proj_name == None or args.project != None:
            self.proj_name = args.project
        if self.proj_name == None:
            print "Please enter a '--project' to run in."
            sys.exit(1)

        if args.sa != None:
            if isinstance(args.sa,list):
                self.statuses_accepted.extend(args.sa)
            else:
                self.statuses_accepted.append(args.sa)

        self.project = dx.get_project(self.proj_name)
        self.proj_id = self.project.get_id()

        # Resolve exp/replicate parameters
        if args.br != None and args.br != 0:
            if len(args.experiments) != 1:
                print "May only specify a replicate for a single experiment."
                sys.exit(1)
            self.rep = { 'br': args.br, 'tr': args.tr, 'rep_tech': 'rep' + str(args.br) + '_' + str(args.tr) }
        else:
            self.rep = None

        if args.genome:
            self.genome = args.genome


    def find_replicates(self, exp_id, exp, verbose=False):
        '''Returns a list of replicates with input files for this experiment.'''
        if self.rep != None:
            return [ self.rep ]
        #verbose=True

        # Must look through exp and find all replicates!
        if exp != self.exp or self.full_mapping == None:
            self.full_mapping = encd.get_full_mapping(exp_id,exp)
        replicates = encd.get_reps(exp_id, exp=exp, full_mapping=self.full_mapping)
        if verbose:
            print "Replicates from encoded:"
            print json.dumps(replicates,indent=4,sort_keys=True)
        for rep in replicates:
            if self.genome == None:
                if rep['organism'] in dx.GENOME_DEFAULTS:
                    self.genome = dx.GENOME_DEFAULTS[rep['organism']]
                else:
                    print "Organism %s not currently supported" % rep['organism']
                    sys.exit(1)
        #replicates = []
        #for (br,tr) in self.full_mapping.keys():
        #    replicates.append( { 'br': br, 'tr': tr,'rep_tech': 'rep' + str(br) + '_' + str(tr) } )
        #
        #    mapping = encd.get_replicate_mapping(exp_id,br,tr,self.full_mapping)
        #    if self.genome == None:
        #        if mapping['organism'] in encd.GENOME_DEFAULTS:
        #            self.genome = dx.GENOME_DEFAULTS[mapping['organism']]
        #        else:
        #            print "Organism %s not currently supported" % mapping['organism']
        #            sys.exit(1)
        #    elif self.genome != dx.GENOME_DEFAULTS[mapping['organism']]:
        #        print "Mixing genomes in one assembly run not supported %s and %s" % \
        #                                            (self.genome, dx.GENOME_DEFAULTS[mapping['organism']])
        #        sys.exit(1)


        if verbose:
            print "Replicates:"
            print json.dumps(replicates,indent=4)
        return replicates

    def find_encoded_files(self, exp, exp_files,replicates,verbose=False):
        '''Returns list of enc file_objects for all files available on encoded.'''
        enc_files = []
        enc_file_names = []
        enc_file_md5sums = []
        #verbose = True
        rep_techs = []
        for rep in replicates:
            rep_techs.append(rep['rep_tech'])

        # Input file should match on format and have format
        input_files = encd.files_to_map(self.exp)
        if verbose or len(self.statuses_accepted) > len(self.FILE_STATUSES_ACCEPTED):
            print ". Accepted file statuses:"
            print self.statuses_accepted
        for f_obj in input_files:
            if f_obj.get('status') not in self.statuses_accepted:
                continue
            rep_tech = 'rep' + str(f_obj['replicate']['biological_replicate_number']) + '_' + str(f_obj['replicate']['technical_replicate_number'])
            if rep_tech not in rep_techs:
                continue
            file_md5 = f_obj['md5sum']
            if file_md5 in enc_file_md5sums:
                continue
            if f_obj.get('file_format') == 'fastq':
                #f_obj['dx_file_name'] = f_obj['accession'] + ".fastq.gz'
                enc_files.append( f_obj )
                enc_file_md5sums.append(file_md5)

        # Result file must match their glob!
        if not self.inputs_only:
            result_files = encd.get_exp_files(exp)
            for obj_type in exp_files['results'].keys():
                for f_obj in result_files:
                    if verbose:
                        print ". Considering %s output_type '%s''" % (f_obj['submitted_file_name'],f_obj['output_type'])
                    bio_rep = f_obj.get('biological_replicates',[])
                    tech_rep = f_obj.get('technological_replicates',[])
                    if len(bio_rep) == 1 and len(tech_rep) == 1:
                        rep_tech = 'rep' + str(bio_rep[0]) + '_' + str(tech_rep[0])
                        if rep_tech not in rep_techs:
                            if verbose:
                                print "~ Excluding because rep_tech is " + rep_tech
                            continue
                    rep = f_obj.get('replicate')
                    if rep != None:
                        rep_tech = 'rep' + str(rep['biological_replicate_number']) + '_' + str(rep.get('technical_replicate_number',1))
                        if rep_tech not in rep_techs:
                            if verbose:
                                print "~ Excluding because 'replicate' rep_tech is " + rep_tech
                            continue
                    if self.genome is not None and 'assembly' in f_obj:
                        if self.genome != f_obj['assembly']:
                            if verbose:
                                print "~ Skipping file %s as its %s != %s." % (f_obj.get('accession'),f_obj['assembly'],self.genome)
                            continue
                    if obj_type != f_obj['output_type']:
                            if verbose:
                                print "~ Skipping file %s as its opbject type '%s'' != '%s'." % (f_obj.get('accession'),f_obj['output_type'],obj_type)
                            continue
                    if f_obj.get('award') != encd.DEFAULT_DCC_AWARD:
                        if verbose:
                            print "~ Skipping file %s as it is not from DCC pipeline." % f_obj.get('accession')
                        continue
                    file_path = f_obj['submitted_file_name'] # Can check submitted_file_name because these should
                    if file_path in enc_file_names:          # have been submitted from dx pipeline
                        #if verbose:
                        #    print "~ Already found."
                        continue
                    for glob in exp_files['results'][obj_type]:
                        if fnmatch.fnmatch(file_path, glob):
                            enc_files.append( f_obj )
                            enc_file_names.append(file_path)
                            break
                        elif verbose:
                            print "~ Not glob: " + glob

        if verbose:
            print ". Encoded files:"
            #print json.dumps(enc_files,indent=4)
            for f_obj in enc_files:
                print ". %s %s" % (f_obj['href'],f_obj['submitted_file_name'])
        return enc_files


    def subtract_files_already_in_dx(self, enc_files, exp_id, exp_folder,verbose=False):
        '''Subtracts from an enc file_objs list any files that are already in the dx exp_folder.'''
        needed_files = []
        for f_obj in enc_files:
            # Inputs are handled one way:
            if f_obj.get('file_format') == 'fastq':
                dx_file_name = os.path.basename(f_obj['href'])
                #file_path = exp_folder + dx_file_name
                file_path = dx_file_name
                # Input files may be at exp folder level!
                # Actually input files can be found anywhere in the project, so use recurse
                #fid = dx.find_file(exp_folder + dx_file_name,self.proj_id,recurse=True)
                fid = dx.find_file(dx_file_name,self.proj_id,recurse=True)
                if fid == None:
                    f_obj['dx_file_name'] = dx_file_name
                    br = f_obj['replicate']['biological_replicate_number']
                    tr = f_obj['replicate']['technical_replicate_number']
                    rep_tech = "rep%d_%d" % (br,tr)
                    f_obj['dx_folder'] = exp_folder + rep_tech + '/'
                    needed_files.append(f_obj)
            else:
                file_path = f_obj['submitted_file_name']
                dx_file_name = os.path.basename(file_path)
                rep_tech = None
                folders = file_path.split('/')
                while folders:
                    folder = folders.pop(0)
                    if folder == exp_id:
                        if folders:
                            folder = folders.pop(0)
                            if folder.startswith('rep') and folder.find('_') != -1:
                                (br,tr) = folder[3:].split('_')
                                if int(br) > 0 and int(tr) > 0:
                                    rep_tech = folder
                        break
                # Now check that encoded.rep_tech hasn't changed.
                enc_rep_tech = rep_tech
                trs = f_obj.get('technical_replicates',[])
                if len(trs) == 1:
                    enc_rep_tech = 'rep' + trs[0]

                # Now look if the file is in dx
                dx_folder = exp_folder
                if rep_tech != None:
                    dx_folder = exp_folder + rep_tech + '/'
                fid = dx.find_file(dx_folder + dx_file_name,self.proj_id,recurse=False)
                if fid == None and enc_rep_tech != rep_tech:
                    dx_folder = exp_folder + enc_rep_tech + '/'
                    fid = dx.find_file(dx_folder + dx_file_name,self.proj_id,recurse=False)
                if fid == None:
                    f_obj['dx_file_name'] = dx_file_name
                    f_obj['dx_folder']    = dx_folder
                    needed_files.append(f_obj)
                elif enc_rep_tech != rep_tech:
                    print "WARNING: found %s in old DX folder: %s" % (dx_file_name,dx_folder)
        if verbose:
            print "Encoded files:"
            print json.dumps(needed_files,indent=4)
        return needed_files


    def prepare_files_to_fetch_json(self, needed_files,verbose=False):
        '''Prepares a json string for requesting files to fetch from encoded to dnanexus.'''
        f2f_files = []
        for f_obj in needed_files:
            f2f_obj = {}     # { "accession": ,"dx_folder": ,"dx_file_name": }
            f2f_obj['accession'] = f_obj['accession']
            f2f_obj['dx_folder'] = f_obj['dx_folder']
            f2f_obj['dx_file_name'] = f_obj['dx_file_name']
            (enc_file_name, bucket_url) = encd.get_bucket(f_obj)
            f2f_obj['enc_file_name'] = enc_file_name
            f2f_obj['bucket_url'] = bucket_url
            f2f_files.append(f2f_obj)
        f2f_json = json.dumps(f2f_files)
        if verbose:
            print "Files to fetch from encoded to dnanexus files json:"
            print f2f_json
        return f2f_json


    def fetch_to_dx(self,exp_id,dx_folder,needed_files,test=True):
        '''Runs fetch-to-dx app to fetch of all files for a given experiment from encoded to dnanexus.'''
        needed_count = len(needed_files)
        files_to_fetch = self.prepare_files_to_fetch_json(needed_files,verbose=False)
        assert (files_to_fetch != None)

        if test:
            print "  - Test fetch %d files from encoded:%s to dnanexus:%s" % (needed_count,exp_id,dx_folder)
            return 0 # Returns the number of files NOT successfully fetched
        else:
            applet = dx.find_applet_by_name('fetch-to-dx', self.proj_id )
            job = applet.run({
                "exp_acc": exp_id,
                "files_to_fetch": files_to_fetch,
                "key": self.server_key,
                "skipvalidate": True,
                "debug": True
                },
                folder=dx_folder,name="Fetch files for "+exp_id)
            print "  - Fetching %d files from encoded:%s to dnanexus:%s  (job:%s)" % \
                    (needed_count,exp_id,dx_folder,job.id)
            sys.stdout.flush() # Slow running job should flush to piped log
            try:
                job.wait_on_done(interval=3)
            except Exception as e:
                print "  " + e.message
                return needed_count

            job_dict = job.describe()
            #error = job_dict['output'].get('error', None)
            if job_dict["state"] == "done":
                fetched_count = job_dict['output'].get('fetched_count', 0)
                return needed_count - fetched_count

        return needed_count # Returns the number of files NOT successfully fetched


    def launch(self,exp_id,exp_type,replicates,genome,annotation,test=True,verbose=False):
        '''
        Spawns the appropriate launcher, not waiting around for the results.
        Returns pid, 0 for noop and -1 for error.
        '''
        # NOT EXPECTED TO OVERRIDE

        # look up the launcher launcher
        if exp_type not in self.LAUNCHERS:
            print "ERROR: No launcher defined for experiment of type " + exp_type
            return -1
        cmd = [ self.LAUNCHERS[exp_type] ]
        cmd.append('-e')
        cmd.append(exp_id)

        # determine arguments to launcher
        # do anything about genome?  Defaults to hg19 or mm10 so not needed yet
        if genome not in ['hg19','mm10']:
            cmd.append('--genome')
            cmd.append(genome)

        # do anything about annotation (lrna, rampage)?  Defaults to v19 or M4 which is okay for now
        if annotation and annotation not in ['v19','M4']:
            cmd.append('--annotation')
            cmd.append(annotation)

        # TODO: Any additional arguments that a launcher might need.

        # Always run, because assemble --test will not actually spawn the command.
        cmd.append('--run')

        if verbose:
            print "Launch command:"
            print json.dumps(cmd,indent=4)

        echo_cmd = ['echo','"']
        echo_cmd.extend(cmd)
        echo_cmd.append('"')
        sys.stdout.flush() # Slow running job should flush to piped log
        if test: # Wait for results and print them
            # NOTE: It would be nice to wait for results of test command, but the command is going to fail since
            #       files have not actually been assembled on a test!
            # Instead, the best we can do is just print the command that would run:
            print "  - Would ignite launcher as: "
            subprocess.call(echo_cmd)
            return 0
        else: # spawn the demon child
            print "  - Ignite launcher as: "
            subprocess.call(echo_cmd)
            subprocess.call(['mkdir','-p','logs/launch/'])
            log_file = 'logs/launch/' + exp_id + '.log'
            with open(log_file,"a") as out:  # append log
                return subprocess.Popen(cmd,stdout=out,stderr=subprocess.STDOUT).pid

    def run(self):
        '''Runs assemble from start to finish using command line arguments.'''
        # NOT EXPECTED TO OVERRIDE

        #try:
        args = self.get_args()
        self.load_variables(args)
        #except Exception as e:
        #    print 'Caught: %s %s' % (e.status_code, e.reason)
        #    #raise

        print "== Running in project [%s] and will copy from the [%s] server ==" % \
                                                        (self.proj_name,self.server_key)

        exp_count = 0
        skipped = 0
        total_copied = 0
        total_failed = 0
        total_launched = 0
        for exp_id in args.experiments:
            sys.stdout.flush() # Slow running job should flush to piped log
            exp_count += 1
            # 1) Lookup experiment type from encoded, based on accession
            self.exp_id = exp_id
            self.exp = encd.get_exp(self.exp_id,must_find=False)
            if self.exp == None or self.exp["status"] == "error":
                print "ERROR: Unable to locate experiment %s in encoded" % exp_id
                skipped += 1
                continue
            if 'internal_status' in self.exp and self.exp['internal_status'] in encd.INTERNAL_STATUS_BLOCKS:
                print "ERROR: Experiment %s with internal_status of '%s' cannot be assembled." % \
                                                                        (self.exp_id,self.exp['internal_status'])
                skipped += 1
                continue
            #print json.dumps(self.exp['files'],indent=4)
            #sys.exit(1)
            self.exp_type = encd.get_exp_type(exp_id,self.exp,self.EXPERIMENT_TYPES_SUPPORTED)
            if self.exp_type == None:
                skipped += 1
                continue

            # Find requested replicate or all replicates for experiment
            self.replicates = self.find_replicates(self.exp_id,self.exp)
            if self.replicates == None or len(self.replicates) == 0:
                print "Warning: No replicates found in encoded for %s." % self.exp_id
                skipped += 1
                continue

            # Ready to announce:
            descr = "ENC2"
            if self.exp.get('award') != None and self.exp['award'].get("rfa") == "ENCODE3":
                descr = "ENC3"
            if self.genome != None:
                descr += "  " +self.genome
            else:
                descr += "  hg19"
            if args.annotation != None:
                descr += " " + args.annotation + " "
            else:
                descr += " v19 "
            if 'shRNA' in self.exp.get('assay_term_name'):
                descr += ' shRNA'
            else:
                descr += '      '
            for rep in self.replicates:
                if descr[-1] != ' ':
                    descr += ','
                descr += str(rep['br'])
            print "Handling %s  %s  %s experiment..." % (exp_id, descr, self.exp_type)

            # File list from encoded
            expected_files = self.PIPELINE_FILES[self.exp_type]
            available_files = self.find_encoded_files(self.exp, expected_files,self.replicates)
            if len(available_files) == 0:
                print "Warning: No files are available in encoded for %s." % self.exp_id
                skipped += 1
                continue
            print "There are %s files available in encoded for %s." % (len(available_files),self.exp_id)

            # 2) Locate the experiment accession named folder
            # NOTE: genome and annotation may have been entered as args to help organize folders
            # NOTE2: annotation may be relevant to launching it is no longer desired for folders
            if self.genome == "mm10" and args.annotation is None:
                args.annotation = "M4"
            self.umbrella_folder = dx.umbrella_folder(args.folder,self.FOLDER_DEFAULT,self.proj_name, \
                                                                            self.exp_type,"runs/",self.genome)
            self.exp_folder = dx.find_exp_folder(self.project,self.exp_id,self.umbrella_folder)
            if self.exp_folder == None:
                self.exp_folder = self.umbrella_folder + exp_id + '/'
                # create it!
                if not self.test:
                    self.project.new_folder(self.exp_folder,parents=True)
                    print "- Have created folder: %s:%s" % (self.proj_name, self.exp_folder)
                else:
                    print "- Test created folder: %s:%s" % (self.proj_name, self.exp_folder)
            else:
                print "- Will examine folder: %s:%s" % (self.proj_name, self.exp_folder)
            for rep in self.replicates:
                rep_folder = self.exp_folder + rep['rep_tech'] + '/'
                if not dx.project_has_folder(self.project, rep_folder):
                    if not self.test:
                        self.project.new_folder(rep_folder)
                        print "  - Have created rep folder:      " + rep_folder
                    else:
                        print "  - Test created rep folder:      " + rep_folder
                else:
                    print "  - Will examine rep folder:      " + rep_folder

            # minus file list already in dx
            needed_files = self.subtract_files_already_in_dx(available_files, self.exp_id,self.exp_folder)
            if len(needed_files) == 0:
                print "* No files need to be copied to dx for %s" % self.exp_id
                skipped += 1
                continue
            ## short circuit for test
            #needed_files = needed_files[0:1]
            print "- Need to copy %d files to dx for %s" % (len(needed_files),self.exp_id)

            # Now for each needed report it
            already_processed_files = False
            for f_obj in needed_files:
                sys.stdout.flush() # Slow running job should flush to piped log
                if not self.test:
                    print "  - Will try to copy %s to dx %s:%s%s" % \
                            (f_obj['accession'],self.proj_name,f_obj['dx_folder'],f_obj['dx_file_name'])
                else:
                    print "  - Would try to copy %s to dx %s:%s%s" % \
                            (f_obj['accession'],self.proj_name,f_obj['dx_folder'],f_obj['dx_file_name'])
                if not f_obj['dx_file_name'].endswith('.fastq.gz') and not f_obj['dx_file_name'].endswith('.fq.gz'):
                    already_processed_files = True

            failed = 0
            failed = self.fetch_to_dx(self.exp_id,self.exp_folder,needed_files,test=self.test)
            copied = len(needed_files) - failed
            launched = 0

            if failed == 0 and not already_processed_files:
                encd.exp_patch_internal_status(self.exp_id, 'pipeline ready', test=self.test)
            #else:
            #    encd.exp_patch_internal_status(self.exp_id, 'unrunnable', test=self.test)

            # Ignite a launcher here...
            if args.launch:
                pid = self.launch(self.exp_id,self.exp_type,self.replicates,self.genome,args.annotation,test=self.test)
                if pid > 0:
                    #print "- Launcher ignited for %s, pid %d" % (self.exp_id, pid)
                    launched = 1

            if self.test:
                print "- For %s Processed %d file(s), would try to copy %d file(s)" % \
                                                            (self.exp_id, len(needed_files), copied)
            else:
                print "- For %s Processed %d file(s), copied %d, failed %d, launched %d" % \
                                                            (self.exp_id, len(needed_files), copied, failed, launched)

            total_copied += copied
            total_failed += failed
            total_launched += launched

        if exp_count > 1:
            if self.test:
                print "Processed %d experiment(s), skipped %d, would try to copy %d file(s)" % \
                                                                (exp_count, skipped, total_copied)
            else:
                print "Processed %d experiment(s), skipped %d, copied %d file(s), failures %d, launched %d" % \
                                                            (exp_count, skipped, total_copied, total_failed, total_launched)
        if total_failed != 0:
            print "(finished with failures)"
            sys.exit(1)

        print "(finished)"

if __name__ == '__main__':
    '''Run from the command line.'''
    # EXPECTED to crash and burn on base class!
    assemble = Assemble()
    assemble.run()

