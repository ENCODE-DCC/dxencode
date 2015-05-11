#!/usr/bin/env python
# assemble.py 0.0.2

import argparse,os, sys, json
#import urlparse, subprocess, itertools, logging
#from datetime import datetime
#import Lib/fnmatch.py
import fnmatch

import dxpy
import dxencode

class Assemble(object):
    '''
    Assembles all input files in dx for a given accessioned experiment and replicate (or combined replicates).  
    '''

    SERVER_DEFAULT = 'www'
    '''At this time there is no need to use the any but the one true server for assembling.'''

    FOLDER_DEFAULT = '/runs/'
    '''This the default location for creating experiment folders on dnanexus.'''
    
    EXPERIMENT_TYPES_SUPPORTED = [ 'long-rna-seq', 'small-rna-seq', 'rampage', 'dnase-seq' ] #,"dnase" ,"dna-me","chip-seq" ]
    '''This module supports only these experiment (pipeline) types.'''

    # Pipeline files includes inputs and results.  To assemble the files, there is no need to understand step order
    # dependencies.  Even replicate vs. combined is not important information as folder destinations can be discovered
    # from submitted_file_name paths. 
    # For each "output_type" there will be one or more globs to recognize files.  Input files really don't need globs.
    PIPELINE_FILES = {
         "long-rna-seq": {
            "inputs":  { "reads": [ "*.fastqs.gz" ] },
            "results": { "alignments":                 [ "*_star_genome.bam",         "*_tophat.bam"          ],
                         "multi-read signal":          [ "*_star_genome_all.bw",      "*_tophat_all.bw"       ],
                         "unique signal":              [ "*_star_genome_uniq.bw",     "*_tophat_uniq.bw"      ],
                         "multi-read minus signal":    [ "*_star_genome_minusAll.bw", "*_tophat_minusAll.bw"  ],
                         "multi-read plus signal":     [ "*_star_genome_plusAll.bw",  "*_tophat_plusAll.bw"   ],
                         "unique minus signal":        [ "*_star_genome_minusUniq.bw","*_tophat_minusUniq.bw" ],
                         "unique plus signal":         [ "*_star_genome_plusUniq.bw", "*_tophat_plusUniq.bw"  ],
                         "transcriptome alignments":   [ "*_star_anno.bam"          ],
                         "genome quantifications":     [ "*_rsem.genes.results"     ],
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
            "results": { "alignments":                 [ "*_rampage_star_marked.bam" ],
                         "multi-read minus signal":    [ "*_rampage_5p_minusAll.bw"  ],
                         "multi-read plus signal":     [ "*_rampage_5p_plusAll.bw"   ],
                         "unique minus signal":        [ "*_rampage_5p_minusUniq.bw" ],
                         "unique plus signal":         [ "*_rampage_5p_plusUniq.bw"  ],
                         "sites":                      [ "*_rampage_idr.gff"         ],  # TODO: not really "sites"
                         "peaks":                      [ "*_rampage_idr.bb"          ] } # TODO: not really "peaks" ?
            # TODO: "controls"?  Punt at this time
            # TODO: "references" Probably not.  Assume that they are in place.
        },
        "dnase-seq": { # TODO: Flesh out the DNase results
            "inputs":  { "reads": [ "*.fastqs.gz" ] },
            "results": { "alignments":                 [ "*_bwa.bam"    ] }
        }
    }
    
    GENOMES_SUPPORTED = ['hg19', 'mm10']
    GENOME_DEFAULT = 'hg19'
    ''' This the default Genome that long RNA-seq experiments are mapped to.'''

    def __init__(self):
        '''
        Assemble expects to be the only class for assembling experiments on for pipeline types.
        '''
        self.args = {} # run time arguments
        self.server_key = self.SERVER_DEFAULT
        self.proj_name = None
        self.project = None
        self.proj_id = None
        self.exp = {}  # Will hold the encoded exp json
        self.psv = {} # will hold pipeline specific variables.
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
                                                        dxencode.env_get_current_project() + "')",
                        required=False)

        # Don't anticpate dealing with reference file
        #ap.add_argument('--refLoc',
        #                help="The location to find reference files (default: '" + \
        #                    dxencode.REF_PROJECT_DEFAULT + ":" + dxencode.REF_FOLDER_DEFAULT + "')",
        #                default=dxencode.REF_FOLDER_DEFAULT,
        #                required=False)

        ap.add_argument('-g','--genome',
                        help='Optionally enter genome to help organize folders',
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

        ap.add_argument('--test',
                        help='Test run only, do not launch anything.',
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
        if self.server_key != "test":
            self.acc_prefix = "ENCFF"
        self.proj_name = dxencode.env_get_current_project()
        if self.proj_name == None or args.project != None:
            self.proj_name = args.project
        if self.proj_name == None:
            print "Please enter a '--project' to run in."
            sys.exit(1)

        self.project = dxencode.get_project(self.proj_name)
        self.proj_id = self.project.get_id()

        # Resolve exp/replicate parameters
        if args.br != None and args.br != 0:
            if len(args.experiments) != 1:
                print "May only specify a replicate for a single experiment." 
                sys.exit(1)
            self.rep = { 'br': args.br, 'tr': args.tr, 'rep_tech': 'rep' + str(args.br) + '_' + str(args.tr) }
        else:
            self.rep = None
        
        
    def find_replicates(self, exp_id, exp, verbose=False):
        '''Returns a list of replicates with input files for this experiment.'''
        if self.rep != None:
            return [ self.rep ]
            
        # Must look through exp and find all replicates!
        self.full_mapping = dxencode.get_full_mapping(exp_id,exp,key=self.server_key)
        replicates = []
        for (br,tr) in self.full_mapping.keys():
            replicates.append( { 'br': br, 'tr': tr,'rep_tech': 'rep' + str(br) + '_' + str(tr) } )
            
        if verbose:
            print "Replicates:"
            print json.dumps(replicates,indent=4)
        return replicates

    def find_encoded_files(self, exp, exp_files,verbose=False):
        '''Returns list of enc file_objects for all files available on encoded.'''
        enc_files = []
        enc_file_names = []
        
        # Input file should match on format and have format
        input_files = dxencode.files_to_map(self.exp)
        for f_obj in input_files:
            if f_obj.get('status') not in ["released","uploaded"]:
                continue
            file_path = f_obj['submitted_file_name']
            if file_path in enc_file_names:
                continue
            if f_obj.get('file_format') == 'fastq':
                #f_obj['dx_file_name'] = f_obj['accession'] + ".fastq.gz'
                enc_files.append( f_obj ) 
                enc_file_names.append(file_path) 

        # Result file must match their glob!
        if not self.inputs_only:
            result_files = dxencode.get_enc_exp_files(exp,key=self.server_key)
            for obj_type in exp_files['results'].keys():
                for f_obj in result_files:
                    if obj_type == f_obj['output_type']:
                        file_path = f_obj['submitted_file_name']
                        if file_path in enc_file_names:
                            continue
                        for glob in exp_files['results'][obj_type]:
                            if fnmatch.fnmatch(file_path, glob):
                                enc_files.append( f_obj )
                                enc_file_names.append(file_path) 

        if verbose:
            print "Encoded files:"
            #print json.dumps(enc_files,indent=4)
            for f_obj in enc_files:
                print "%s %s" % (f_obj['href'],f_obj['submitted_file_name'])
        return enc_files
        
        
    def subtract_files_already_in_dx(self, enc_files, exp_id, exp_folder,verbose=False):
        '''Subtracts from an enc file_objs list any files that are already in the dx exp_folder.'''
        needed_files = []
        for f_obj in enc_files:
            # Inputs are handled one way:
            if f_obj.get('file_format') == 'fastq':
                dx_file_name = os.path.basename(f_obj['href'])
                file_path = exp_folder + dx_file_name
                # Input files may be at exp folder level!
                # Actually input files can be found anywhere in the project, so use recurse
                fid = dxencode.find_file(exp_folder + dx_file_name,self.proj_id,recurse=True)
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
                # Now look if the file is in dx
                dx_folder = exp_folder
                if rep_tech != None:
                    dx_folder += rep_tech + '/'
                fid = dxencode.find_file(dx_folder + dx_file_name,self.proj_id,recurse=False)
                if fid == None:
                    f_obj['dx_file_name'] = dx_file_name
                    f_obj['dx_folder']    = dx_folder
                    needed_files.append(f_obj)
        if verbose:
            print "Encoded files:"
            print json.dumps(needed_files,indent=4)
        return needed_files


    def prepare_files_to_fetch_json(self, needed_files,verbose=False):
        '''Prepares a json string for requesting files to fetch from encoded to dnanexus.'''
        f2f_files = []
        AUTHID, AUTHPW, SERVER = dxencode.processkey(self.server_key)
        for f_obj in needed_files:
            f2f_obj = {}     # { "accession": ,"dx_folder": ,"dx_file_name": }
            f2f_obj['accession'] = f_obj['accession']
            f2f_obj['dx_folder'] = f_obj['dx_folder']
            f2f_obj['dx_file_name'] = f_obj['dx_file_name']
            (enc_file_name, bucket_url) = dxencode.get_bucket(SERVER, AUTHID, AUTHPW, f_obj)
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
            applet = dxencode.find_applet_by_name('fetch-to-dx', self.proj_id )
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


    def run(self):
        '''Runs launch from start to finish using command line arguments.'''
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
        for exp_id in args.experiments:
            sys.stdout.flush() # Slow running job should flush to piped log
            exp_count += 1
            # 1) Lookup experiment type from encoded, based on accession
            self.exp_id = exp_id
            self.exp = dxencode.get_exp(self.exp_id,must_find=False,key=self.server_key)
            if self.exp == None or self.exp["status"] == "error":
                print "Unable to locate experiment %s in encoded" % exp_id
                skipped += 1
                continue
            #print json.dumps(self.exp['files'],indent=4)
            #sys.exit(1)
            self.exp_type = dxencode.get_exp_type(exp_id,self.exp,self.EXPERIMENT_TYPES_SUPPORTED)
            if self.exp_type == None:
                skipped += 1
                continue
            expected_files = self.PIPELINE_FILES[self.exp_type]
            print "Handling %s as '%s' experiment..." % (exp_id,self.exp_type)
                
            # Find requested replicate or all replicates for experiment
            self.replicates = self.find_replicates(self.exp_id,self.exp)
            if self.replicates == None or len(self.replicates) == 0:
                print "Warning: No replicates found in encoded for %s." % self.exp_id
                skipped += 1
                continue

            # File list from encoded
            available_files = self.find_encoded_files(self.exp, expected_files)
            if len(available_files) == 0:
                print "Warning: No files are available in encoded for %s." % self.exp_id
                skipped += 1
                continue
            print "There are %s files available in encoded for %s." % (len(available_files),self.exp_id)
            
            # 2) Locate the experiment accession named folder
            # NOTE: genome and annotation may have been entered as args to help organize folders
            self.umbrella_folder = dxencode.umbrella_folder(args.folder,self.FOLDER_DEFAULT,self.exp_type, \
                                                                                            args.genome,args.annotation)
            self.exp_folder = dxencode.find_exp_folder(self.project,self.exp_id,self.umbrella_folder)
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
                if not dxencode.project_has_folder(self.project, rep_folder):
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
            for f_obj in needed_files:
                sys.stdout.flush() # Slow running job should flush to piped log
                if not self.test: 
                    print "  - Will try to copy %s to dx %s:%s%s" % \
                            (f_obj['accession'],self.proj_name,f_obj['dx_folder'],f_obj['dx_file_name'])
                else:
                    print "  - Would try to copy %s to dx %s:%s%s" % \
                            (f_obj['accession'],self.proj_name,f_obj['dx_folder'],f_obj['dx_file_name'])

            failed = 0
            failed = self.fetch_to_dx(self.exp_id,self.exp_folder,needed_files,test=self.test)
            copied = len(needed_files) - failed
            
            if self.test:
                print "- For %s Processed %d file(s), would try to copy %d file(s)" % \
                                                            (self.exp_id, len(needed_files), copied)
            else:
                print "- For %s Processed %d file(s), copied %d, failed %d" % \
                                                            (self.exp_id, len(needed_files), copied, failed)
            total_copied += copied
            total_failed += failed
        
        if exp_count > 1:
            if self.test:
                print "Processed %d experiment(s), skipped %d, would try to copy %d file(s)" % \
                                                                (exp_count, skipped, total_copied)
            else:    
                print "Processed %d experiment(s), skipped %d, copied %d file(s) failures %d" % \
                                                                (exp_count, skipped, total_copied, total_failed)
        if total_failed != 0:
            print "(finished with failures)"
            sys.exit(1)    

        print "(finished)"

if __name__ == '__main__':
    '''Run from the command line.'''
    # EXPECTED to crash and burn on base class!
    assemble = Assemble()
    assemble.run()

