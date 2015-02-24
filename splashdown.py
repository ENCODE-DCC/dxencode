#!/usr/bin/env python
# splashdown.py 0.0.1
#
# Initial starting point accessonator.py in tf_chipseq.py and lrnaSplashdown.py
#
# Splashdown is meant to run outside of dnanexus and to examine experiment directories to
# find results tp post to encoded.
#
# 1) Lookup experiment type from encoded, based on accession
# 2) Locate the experiment accession named folder
# 3) Given the experiment type, determine the expected results
# 4) Given expected results locate any files (by glob) that should be posted for
#    a) each single replicate (in replicate sub-folders named as reN_N/
#    b) combined replicates in the experiment folder itself
# 5) For each file that should be posted, determine if the file needs to be posted (not already).
# 6) For each file that needs to be posted:
#    a) discover all necessary dx information needed for posting.
#    b) gather any other information necessary from dx and encoded. (notice hand-waving)
#    c) Post file and update encoded database.
#    d) Update dnanexus file with file accession tag.
# 7) Either exit or advance to the next experiment folder
 # NOTE: any job output name='metadata' class=string (in json format!) will be added to each postable file of that job

import argparse,os, sys
import json, urlparse, subprocess, itertools, logging
#import requests, re, shlex, time
#from datetime import datetime

import dxpy
import dxencode

class Splashdown(object):
    '''
    Splashdown module posts from dnanexus to ENCODEd,  all files available and necessry for
    a given experiment .
    '''

    SERVER_DEFAULT = 'test'
    '''This the default server to post files to.'''

    RESULT_FOLDER_DEFAULT = "/"
    '''Where to start the search for experiment folders.'''

    EXPERIMENT_TYPES_SUPPORTED = [ 'long-rna-seq', 'small-rna-seq', 'rampage' ] #,"dna-me","chip-seq" ]
    '''This module supports only these experiment (pipeline) types.'''

    SKIP_VALIDATE = []
    '''Some output_types cannot currently be validated, theoretically'''

    #Pipeline specifications include order of steps, steps per replicate, combined steps and
    #within steps, the output_type: file_glob that define expected results.
    PIPELINE_SPECS = {
         "long-rna-seq": {
            "step-order": [ "align-tophat","signals-top-se","signals-top-pe",
                            "align-star","signals-star-se","signals-star-pe","quant-rsem"],
            "replicate":  {
                "align-tophat":    { "alignments":                "*_tophat.bam"                 },
                "signals-top-se":  { "multi-read signal":         "*_tophat_all.bw",
                                     "unique signal":             "*_tophat_uniq.bw"             },
                "signals-top-pe":  { "multi-read minus signal":   "*_tophat_minusAll.bw",
                                     "multi-read plus signal":    "*_tophat_plusAll.bw",
                                     "unique minus signal":       "*_tophat_minusUniq.bw",
                                     "unique plus signal":        "*_tophat_plusUniq.bw"         },
                "align-star":      { "alignments":                "*_star_genome.bam",
                                     "transcriptome alignments":  "*_star_anno.bam"              },
                "signals-star-se": { "multi-read signal":         "*_star_genome_all.bw",
                                     "unique signal":             "*_star_genome_uniq.bw"        },
                "signals-star-pe": { "multi-read minus signal":   "*_star_genome_minusAll.bw",
                                     "multi-read plus signal":    "*_star_genome_plusAll.bw",
                                     "unique minus signal":       "*_star_genome_minusUniq.bw",
                                     "unique plus signal":        "*_star_genome_plusUniq.bw"    },
                "quant-rsem":      { "transcript quantifications":"*_rsem.isoforms.results",
                                     "genome quantifications":    "*_rsem.genes.results"         }  },
            "combined":   {}
        },
        "small-rna-seq": {
            "step-order": [ "align","signals"],
            "replicate":  {
                "align":           { "alignments":                "*_star_genome.bam"            },
                "signals":         { "multi-read plus signal":    "*_small_plusAll.bw",
                                     "multi-read minus signal":   "*_small_minusAll.bw",
                                     "unique plus signal":        "*_small_plusUniq.bw",
                                     "unique minus signal":       "*_small_minusUniq.bw"         }  },
            "combined":   {}
        },
        "rampage": {
            "step-order": [ "align","signals","peaks","idr"],
            "replicate":  {
                "align":           { "alignments":                "*_rampage_star_marked.bam" },
                "signals":         { "multi-read plus signal":    "*_rampage_5p_plusAll.bw",
                                     "multi-read minus signal":   "*_rampage_5p_minusAll.bw",
                                     "unique plus signal":        "*_rampage_5p_plusUniq.bw",
                                     "unique minus signal":       "*_rampage_5p_minusUniq.bw" },
                "peaks":           { "peaks":                     "*_rampage_peaks.bb",
                                     "sites":                     "*_rampage_peaks.gff"       } },
            "combined":   {
                "idr":             { "sites":                     "*_rampage_idr.gff",    # TODO: not really "sites"
                                     "peaks":                     "*_rampage_idr.bb" }  } # TODO: not really "peaks" ?
        }
    }

    ASSEMBLIES_SUPPORTED = { "hg19": "hg19", "hg38": "GRCh37", "mm10": "mm10" }
    '''This module supports only these assemblies.'''

    ANNOTATIONS_SUPPORTED = [ 'V19', 'M2', 'M3', 'M4' ]
    '''This module supports only these annotations.'''

    FORMATS_SUPPORTED = ["bam", "bed", "bedLogR", "bed_bedLogR", "bedMethyl", "bed_bedMethyl",
                         "bigBed", "bigWig", "broadPeak", "bed_broadPeak", "fasta", "fastq",
                         "gtf", "idat", "narrowPeak", "bed_narrowPeak", "rcc", "CEL", "tsv", "csv" ]
    EXTENSION_TO_FORMAT = { "bb":"bigBed", "bw":"bigWig",
                            "fa":"fasta","fq":"fastq","results":"tsv",
                            "gff": "gtf" }
    '''List of supported formats, and means of recognizing with file extensions.'''

    PRIMARY_INPUT_EXTENSION = [ "fastq","fq"]
    '''List of file extensions used to recognize primary inputs to parse accessions.'''


    def __init__(self):
        '''
        Splashdown expects one or more experiment ids as arguments and will find, document
        and post files in the associated directory.
        '''
        self.args = {} # run time arguments
        self.server_key = 'test'
        self.acc_prefix = "TSTFF"
        self.proj_name = None
        self.project = None
        self.proj_id = None
        self.exp = {}  # Will hold the encoded exp json
        self.exp_id = None
        self.exp_type = {}  # Will hold the experiment's assay_type, normalized to known tokens.
        self.genome = None  # TODO: need way to determine genome before any posts occur!
        self.annotation = None  # TODO: if appropriate, need way to determine annotation
        self.pipeline = None # pipeline definitions (filled in when experiment type is known)
        self.replicates = None # lost replicate folders currently found beneath experiment folder
        self.test = True # assume Test until told otherwise
        dxencode.logger = logging.getLogger(__name__) # I need this to avoid some errors
        dxencode.logger.addHandler(logging.StreamHandler()) #logging.NullHandler)


    def get_args(self):
        '''Parse the input arguments.'''
        ### PIPELINE SPECIFIC
        ap = argparse.ArgumentParser(description="Handles splashdown of launched pipeline runs " +
                    "for supported experiment types. Can be run repeatedly and will only try to " +
                    "post result files that have not been previously posted. All results " +
                    "are expected to be in folder /<resultsLoc>/<experiment> and any replicate " +
                    "sub-folders named as " +
                    "<experiment>/rep<biological-replicate>_<technical-replicate>.")
        ### PIPELINE SPECIFIC

        ap.add_argument('-e', '--experiments',
                        help='One or more ENCODED experiment accessions',
                        nargs='+',
                        required=True)

        ap.add_argument('--project',
                        help="Project to run analysis in (default: '" + \
                                                         dxencode.env_get_current_project() + "')",
                        required=False)

        ap.add_argument('--results_folder',
                        help="The location to search for experiment folders (default: " + \
                                                "'<project>:" + self.RESULT_FOLDER_DEFAULT + "')",
                        default=self.RESULT_FOLDER_DEFAULT,
                        required=False)

        ap.add_argument('--server',
                        help="Server to post files to (default: '" + self.SERVER_DEFAULT + "')",
                        default=self.SERVER_DEFAULT,
                        required=False)

        ap.add_argument('--test',
                        help='Test run only, do not launch anything.',
                        action='store_true',
                        required=False)

        ap.add_argument('--verbose',
                        help='More debugging output.',
                        action='store_true',
                        required=False)

        ap.add_argument('--force_annotation',
                        help='If annotation cannot be found, use this string',
                        required=False)

        return ap.parse_args()

    def get_exp_type(self,exp_id,exp=None):
        '''Looks up encoded experiment's assay_type, normalized to known supported tokens.'''
        self.exp_id = exp_id
        if exp == None and self.exp == None:
            self.exp = get_exp(exp_id)
        self.exp_type = dxencode.get_assay_type(self.exp_id,self.exp)

        if self.exp_type not in self.EXPERIMENT_TYPES_SUPPORTED:
            print "Experiment %s has unsupported assay type of '%s'" % \
                                                            (exp_id,self.exp["assay_term_name"])
            return None
        return self.exp_type


    def find_replicate_folders(self,exp_folder,verbose=False):
        '''Returns a sorted list of replicate folders which are sub-folders of an exp folder.'''
        # normalize
        try:
            sub_folders = self.project.list_folder(exp_folder)['folders']
        except:
            if verbose:
                print "No subfolders found for %s" % exp_folder
            return []

        replicates = []
        for path in sorted(sub_folders):
            folder = path.split('/')[-1]
            if folder.startswith('rep'):
                if len(folder[3:].split('_')) == 2:
                    replicates.append( folder )
        if verbose:
            print "Replicate folders:"
            print json.dumps(replicates,indent=4)
        return replicates


    def pipeline_specification(self,args,exp_type,exp_folder,verbose=False):
        '''Sets the pipeline specification object for this experiment type.'''

        # Start with dict containing common variables
        #self.expected = copy.deepcopy(self.PIPELINE_SPECS[exp_type])

        pipeline_specs = self.PIPELINE_SPECS.get(exp_type)
        self.genome = None  # TODO: need way to determine genome before any posts occur!
        self.annotation = None  # TODO: if appropriate, need way to determine annotation


        if verbose:
            print "Pipeline specification:"
            print json.dumps(pipeline_specs,indent=4)
        return pipeline_specs


    def file_format(self,file_name):
        '''Try to determine file format from file name extension.'''
        ext = file_name.split(".")[-1]
        if ext == "gz" or ext == "tgz":
            ext = file_name.split(".")[-2]
        if ext in self.EXTENSION_TO_FORMAT.keys():
            ext = self.EXTENSION_TO_FORMAT[ext]
        if ext in self.FORMATS_SUPPORTED:
            return ext
        return None


    def find_genome_annotation(self,file_dict):
        '''Try to determine genome from input file properties.'''
        # TODO: currently done in derived_from which is only run on needed files
        #       much change to do this on expected files.
        properties = file_dict["properties"]
        msg = ""
        if self.genome == None and "genome" in properties:
            genome = properties["genome"]
            if genome in self.ASSEMBLIES_SUPPORTED.keys():
                self.genome = self.ASSEMBLIES_SUPPORTED[genome]
                msg += " genome[%s]" % self.genome
        if self.annotation == None and "annotation" in properties:
            annotation = properties["annotation"].upper() # BRITTLE: v19 in dx but V19 in encoded
            if annotation in self.ANNOTATIONS_SUPPORTED:
                self.annotation = annotation
                msg += " annotation[%s]" % self.annotation
        if len(msg) > 0:
            print "  - Found " + msg
        return self.genome


    def find_step_files(self,file_globs,result_folder,rep_tech,verbose=False):
        '''Returns tuple list of (type,rep_tech,fid) of ALL files expected for a single step.'''
        step_files = []

        for token in file_globs.keys():
            if self.file_format(file_globs[token]) == None:
                print "Error: file glob %s has unknown file format! Please fix" % file_globs[token]
                sys.exit(1)
            if verbose:
                print "-- Looking for %s" % (result_folder + file_globs[token])
            fid = dxencode.find_file(result_folder + file_globs[token],self.proj_id, recurse=False)
            if fid != None:
                step_files.append( (token,rep_tech,fid) )
            else:
                return []      # Only include files from completed steps!
        return step_files

    def find_expected_files(self,exp_folder,replicates,verbose=False):
        '''Returns tuple list of (type,rep_tech,fid) of files expected to be posted to ENCODE.'''
        expected = []
        # First find replicate step files
        for step in self.pipeline["step-order"]:
            if step not in self.pipeline["replicate"]:
                continue
            for rep_tech in replicates:
                step_files = self.find_step_files(self.pipeline["replicate"][step], \
                                                    exp_folder + rep_tech + '/',rep_tech,verbose)
                if len(step_files) > 0:
                     expected.extend(step_files) # keep them in order!

        # Now add combined step files
        for step in self.pipeline["step-order"]:
            if step not in self.pipeline["combined"]:
                continue
            step_files = self.find_step_files(self.pipeline["combined"][step], \
                                                                    exp_folder,"combined",verbose)
            if len(step_files) > 0:
                 expected.extend(step_files) # keep them in order!

        if verbose:
            print "Expected files:"
            print json.dumps(expected,indent=4)
        return expected


    def find_in_encode(self,fid,verbose=False):
        '''Looks for file in encode with the same submitted file name as the fid has.'''
        file_name = dxencode.file_path_from_fid(fid,projectToo=True).split(':')[-1]
        file_size = dxencode.description_from_fid(fid).get('size')

        # TODO: encoded should allow searching on fid !!!
        authid, authpw, server = dxencode.processkey(self.server_key)
        url = urlparse.urljoin(server,
                'search/?type=file&frame=object&submitted_file_name=%s' % file_name)
        response = dxencode.encoded_get(url,authid,authpw)
        try:
            if verbose:
                print "-- Looking for %s (size:%d) in %s" % (file_name, file_size, server)
            response.raise_for_status()
            if response.json()['@graph']:
                encode_file = response.json()['@graph'][0]
                #logger.info("Found potential duplicate: %s" %(duplicate_item.get('accession')))
                if file_size ==  encode_file.get('file_size'):
                    if verbose:
                        print("%s %s: File sizes match, assuming duplicate." % \
                                                   (str(file_size), encode_file.get('file_size')))
                    return encode_file
                else:
                    if verbose:
                        print("%s %s: File sizes differ, assuming new file." % \
                                                   (str(file_size), encode_file.get('file_size')))
                    return None
            else:
                if verbose:
                    print("No duplicate ... proceeding")
                    print(response.text)
                    print(response.json())
                return None
        except:
            if verbose:
                print('Duplicate accession check failed: %s %s' % (response.status_code, response.reason))
                print(response.text)
            return None

    def find_needed_files(self,files_expected,verbose=False):
        '''Returns the tuple list of files that NEED to be posted to ENCODE.'''
        needed = []
        self.found = {}
        for (out_type, rep_tech, fid) in files_expected:
            # Current strategy is to complete an post before updating the accession field.
            # so existence of accession should mean it is already in encoded.
            fileDict = dxencode.description_from_fid(fid,properties=True)
            acc_key = "accession"
            if self.server_key == 'test':
                acc_key = "test_accession"
            if "properties" in fileDict and acc_key in fileDict["properties"]:
                accession = fileDict["properties"][acc_key]
                if accession.startswith(self.acc_prefix) and len(accession) == 11:
                    continue
            # No accession so try to match in encoded by submit_file_name and size
            f_obj =  self.find_in_encode(fid,verbose)
            if f_obj == None:
                needed.append( (out_type,rep_tech,fid) )
            else:
                self.found[fid] = f_obj

        if verbose:
            print "Needed files:"
            print json.dumps(needed,indent=4)
        return needed


    def find_derived_from(self,fid,job,verbose=False):
        '''Returns list of accessions a file is drived from based upon job inouts.'''
        derived_from = []
        for input in job["input"].values():
            if not type(input) == dict:
                continue # not a file input
            dxlink = input.get("$dnanexus_link")
            if dxlink == None:
                continue
            if not type(dxlink) == dict:
                inp_fid = dxlink
            else:
                inp_fid = dxlink.get("id")
            inp_obj = dxencode.description_from_fid(inp_fid,properties=True)
            if inp_obj != None:
                self.genome = self.find_genome_annotation(inp_obj)
                acc_key = "accession"
                if self.server_key == 'test':
                    acc_key = "test_accession"
                if "properties" in inp_obj and acc_key in inp_obj["properties"]:
                    # older runs don't have
                    accession = inp_obj["properties"][acc_key]
                    if accession.startswith(self.acc_prefix) and len(accession) == 11:
                        derived_from.append(accession)
                else: # if file name is primary input (fastq) and is named as an accession
                    if inp_obj["name"].startswith("ENCFF"): # Not test version 'TSTFF'!
                        parts = inp_obj["name"].split('.')
                        ext = parts[-1]
                        if ext in ["gz","tgz"]:
                            ext = parts[-2]
                        if ext in self.PRIMARY_INPUT_EXTENSION:
                            root = parts[0]
                            for acc in root.split('_'): #usually only one
                                if acc.startswith("ENCFF") and len(acc) == 11:
                                    derived_from.append(acc)
        if verbose:
            print "Derived files: for " + dxencode.file_path_from_fid(fid)
            print json.dumps(derived_from,indent=4)
        return derived_from


    def add_encoded_info(self,obj,rep_tech,fid,verbose=False):
        '''Updates an object with information from encoded database.'''
        obj['lab'] = '/labs/j-michael-cherry/' # self.exp['lab']['@id'] Now hard-coded
        obj['award'] = '/awards/U41HG006992/'  # self.exp['award']['@id']

        # Find replicate info
        if rep_tech.startswith("rep"):
            br_tr = rep_tech[3:]
            (br,tr) = br_tr.split('_')
            full_mapping = dxencode.get_full_mapping(self.exp_id,self.exp)
            mapping = dxencode.get_replicate_mapping(self.exp_id,int(br),int(tr),full_mapping)
            obj['replicate'] = mapping['replicate_id']

        if verbose:
            print "After adding encoded info:"
            print json.dumps(obj,indent=4)
        return obj


    def make_payload_obj(self,out_type,rep_tech,fid,verbose=False):
        '''Returns an object for submitting a file to encode, with all dx info filled in.'''
        payload = {}
        payload['dataset'] = self.exp_id
        payload["output_type"] = out_type
        dx_obj = dxencode.description_from_fid(fid)
        #if verbose:
        #    print "dx_obj:"
        #    print json.dumps(dx_obj,indent=4)
        job = dxencode.job_from_fid(fid)
        #if verbose:
        #    print "job:"
        #    print json.dumps(job,indent=4)
        #applet = dxencode.applet_from_fid(fid)
        payload["file_format"] = self.file_format(dx_obj["name"])
        if payload["file_format"] == None:
            print "Warning: file %s has unknown file format!" % dxencode.file_path_from_fid(fid)
        payload["derived_from"] = self.find_derived_from(fid,job, verbose)
        payload['submitted_file_name'] = dxencode.file_path_from_fid(fid,projectToo=True)
        payload['file_size'] = dx_obj["size"]
        #payload['md5sum'] = calculated_md5 # TODO: Find from file properties???
        if self.genome == None:
            print "Warning: could not determine genome assembly! Add properties to reference files."
            #sys.exit(1)
        else:
            payload['assembly'] = self.genome
        if self.annotation != None:
            payload['genome_annotation'] = self.annotation

        dxFile = dxencode.file_handler_from_fid(fid)
        versions = dxencode.get_sw_from_log(dxFile, '\* (\S+)\s+version:\s+(\S+)')
        notes = dxencode.create_notes(dxFile, versions)
        #notes['qc'] = flagstat_parse(bamqc) ????
        # TODO: find json added to job as a result of returns? ??
        if "metadata" in job["output"]:
           notes.update(json.load(job["output"]["metadata"])) # NOTE: output class=string name='metadata' in json format
        payload['notes'] = json.dumps(notes)

        #print "  - Adding encoded information."
        payload = self.add_encoded_info(payload,rep_tech,fid)

        if verbose:
            print "payload from dx info:"
            print json.dumps(payload,indent=4)
        return payload


    def file_mark_accession(self,fid,accession,test=True):
        '''Adds/replaces accession to a file's properties.'''
        file_handler = dxencode.file_handler_from_fid(fid)
        properties = file_handler.get_properties()
        path = dxencode.file_path_from_fid(fid)
        acc_key = "accession"
        if self.server_key == 'test':
            acc_key = "test_accession"
        if acc_key in properties and properties[acc_key] != accession:
            if properties[acc_key] != accession:
                print "Warning: file %s has accession %s but has been posted as accession %s" % \
                    (path,properties[acc_key],accession)
                #sys.exit(1)
        properties[acc_key] = accession
        if test:
            print "  - Test flag %s with accession='%s'" % (path,accession)
        else:
            file_handler.set_properties(properties)
            print "  - Flagged   %s with accession='%s'" % (path,accession)


    def file_post(self,fid,payload,test=True):
        '''Posts a file to encoded.'''
        path = payload['submitted_file_name'].split(':')[1]
        derived_count = len(payload["derived_from"])
        skip_validate = (payload["output_type"] in self.SKIP_VALIDATE)
        val_msg = ""
        if skip_validate:
            val_msg = " UNVALIDATED"
        if test:
            print "  - Test post %s (%d) to '%s'%s" % (path,derived_count,self.server_key,val_msg)
            if self.server_key == "test":
                return "TSTFF00FAKE"
            return "ENCFF00FAKE"
        else:
            out_folder = self.exp_folder + "posts"
            dxencode.find_or_create_folder(self.project, out_folder)
            applet = dxencode.find_applet_by_name('validate-post', self.proj_id )
            job = applet.run({
                "pipe_file": dxpy.dxlink(fid),
                "file_meta": payload,
                "key": self.server_key,
                "skipvalidate": skip_validate,
                "debug": True
                },
                folder=out_folder)
            print "  - Submitting %s to '%s' (derived_from:%d)%s" % \
                                                     (job.id,self.server_key, derived_count,val_msg)
            sys.stdout.flush() # Slow running job should flush to piped log
            try:
                job.wait_on_done(interval=1)
            except Exception as e:
                print "  " + e.message
                return None

            job_dict = job.describe()
            #error = job_dict['output'].get('error', None)
            if job_dict["state"] == "done":
                accession = job_dict['output'].get('accession', None)
                return accession
            else:
                return None

        return None


    def run(self):
        '''Runs splasdown from start to finish using command line arguments.'''
        args = self.get_args()
        self.test = args.test
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
        print "== Running in project [%s] and will post to the [%s] server ==" % \
                                                        (self.proj_name,self.server_key)

        exp_count = 0
        halted = 0
        total_posted = 0
        for exp_id in args.experiments:
            sys.stdout.flush() # Slow running job should flush to piped log
            # 1) Lookup experiment type from encoded, based on accession
            print "Working on %s..." % exp_id
            self.exp = dxencode.get_exp(exp_id,must_find=False,key=self.server_key)
            if self.exp == None or self.exp["status"] == "error":
                print "Unable to locate experiment %s in encoded" % exp_id
                continue
            self.exp_type = self.get_exp_type(exp_id)
            if self.exp_type == None:
                continue

            # 2) Locate the experiment accession named folder
            self.exp_folder = dxencode.find_exp_folder(self.project,exp_id,args.results_folder,warn=True)
            if self.exp_folder == None:
                continue
            print "- Examining %s:%s for '%s' results..." % \
                                            (self.proj_name, self.exp_folder, self.exp_type)

            # 3) Given the experiment type, determine the expected results
            self.pipeline   = self.pipeline_specification(args,self.exp_type,self.exp_folder)
            self.replicates = self.find_replicate_folders(self.exp_folder, verbose=args.verbose)

            # 4) Given expected results locate any files (by glob) that should be posted for
            #    a) each single replicate (in replicate sub-folders named as reN_N/
            #    b) combined replicates in the experiment folder itself
            files_expected = self.find_expected_files(self.exp_folder,self.replicates, verbose=args.verbose)
            print "- Found %d files that are available to post." % len(files_expected)
            if len(files_expected) == 0:
                continue

            # 5) For each file that should be posted, determine if the file needs to be posted.
            files_to_post = self.find_needed_files(files_expected, verbose=args.verbose)
            print "- Found %d files that need to be posted" % len(files_to_post)
            if len(files_to_post) == 0:
                continue

            # 6) For each file that needs to be posted:
            exp_count += 1
            file_count = 0
            post_count = 0
            for (out_type,rep_tech,fid) in files_to_post:
                sys.stdout.flush() # Slow running job should flush to piped log
                # a) discover all necessary dx information needed for post.
                # b) gather any other information necessary from dx and encoded.
                print "  Handle file %s" % dxencode.file_path_from_fid(fid)
                payload = self.make_payload_obj(out_type,rep_tech,fid, verbose=args.verbose)

                file_count += 1
                # c) Post file and update encoded database.
                accession = self.file_post(fid,payload,args.test)
                if accession == None:
                    print "* HALTING %s - post failure could compromise 'derived_from'" % \
                                                                                    (self.exp_id)
                    halted += 1
                    break

                # d) Update dnanexus file with file accession tag.
                if not args.test:
                    post_count += 1
                self.file_mark_accession(fid,accession,args.test)

                #if file_count >= 5:  # Short circuit for test
                #    break

            print "- For %s Processed %d file(s), posted %s" % \
                                                        (self.exp_id, file_count, post_count)
            total_posted += post_count

        print "Processed %d experiment(s), halted %d, posted %d file(s)" % \
                                                            (exp_count, halted, total_posted)
        if halted == exp_count:
            sys.exit(1)
        print "(finished)"


if __name__ == '__main__':
    '''Run from the command line.'''
    recovery = Splashdown()
    recovery.run()

