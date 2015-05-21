#!/usr/bin/env python2.7
# splashdown.py 0.0.1
#
### TODO:
#   1) The creation of the workflow_run encoded object should be migrated to the launchers.
#   2) Splashdown should update the workflow_run object if it already exists.
#   3) Ideally, splashdown.py could poll encoded for workflow_run objects that have not completed,
#      then determine if the have finished or failed and then begin the process of posting any results that are available.
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
import json, urlparse, subprocess, itertools, logging, time
#import requests, re, shlex, time
from datetime import datetime

import dxpy
import dxencode

class Splashdown(object):
    '''
    Splashdown module posts from dnanexus to ENCODEd,  all files available and necessry for
    a given experiment .
    '''

    SERVER_DEFAULT = 'test'
    '''This the default server to post files to.'''

    FOLDER_DEFAULT = "/"
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
                "signals-star-se": { "multi-read signal":         "*_star_genome_all.bw",
                                     "unique signal":             "*_star_genome_uniq.bw"        },
                "signals-star-pe": { "multi-read minus signal":   "*_star_genome_minusAll.bw",
                                     "multi-read plus signal":    "*_star_genome_plusAll.bw",
                                     "unique minus signal":       "*_star_genome_minusUniq.bw",
                                     "unique plus signal":        "*_star_genome_plusUniq.bw"    },
                "align-star":      { "alignments":                "*_star_genome.bam",
                                     "transcriptome alignments":  "*_star_anno.bam"              },
                "quant-rsem":      { "genome quantifications":    "*_rsem.genes.results",
                                     "transcript quantifications":"*_rsem.isoforms.results"      }  },
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
        self.obj_cache = {} # certain things take time to find or create and are needed multiple times
        logging.basicConfig(format='%(asctime)s  %(levelname)s: %(message)s')
        dxencode.logger = logging.getLogger(__name__ + '.dxe') # I need this to avoid some errors
        dxencode.logger.addHandler(logging.StreamHandler()) #logging.NullHandler)
        print



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

        ap.add_argument('-f','--folder',
                        help="The location to search for experiment folders (default: " + \
                                                "'<project>:" + self.FOLDER_DEFAULT + "')",
                        default=self.FOLDER_DEFAULT,
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

        ap.add_argument('--ignore_properties',
                        help='Ignore DNANexus file properties and try to post, --server == test',
                        action='store_true',
                        required=False)

        return ap.parse_args()


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
            print "  - Found" + msg
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


    def enc_lookup_json(self,path,must_find=False):
        '''Attempts to retrieve an ENCODEd json object.'''
        url = self.server + path + '/?format=json&frame=embedded'
        #print url
        response = dxencode.encoded_get(url, self.authid, self.authpw)
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


    def enc_file_find_find_by_dxid(self,dx_fid):
        '''Finds a encoded 'file' object by dnanexus alias.'''
        file_obj = None
        
        file_alias = 'dnanexus:' + dx_fid
        #if pipe_alias in self.obj_cache:
        #    return self.obj_cache[pipe_alias]
        file_obj = self.enc_lookup_json( 'files/' + file_alias,must_find=False)
        #if v:
        #    self.obj_cache[file_alias] = file_obj
        #    print "  - Found file: '%s'" % file_alias
        return file_obj
       

    def find_needed_files(self,files_expected,verbose=False):
        '''Returns the tuple list of files that NEED to be posted to ENCODE.'''
        needed = []
        self.found = {}
        for (out_type, rep_tech, fid) in files_expected:
            # Current strategy is to complete and post before updating the accession field.
            # so existence of accession should mean it is already in encoded.
            fileDict = dxencode.description_from_fid(fid,properties=True)
            acc_key = dxencode.dx_property_accesion_key(self.server)
            if not self.ignore:
                # check file properties
                if "properties" in fileDict and acc_key in fileDict["properties"]:
                    accession = fileDict["properties"][acc_key]
                    if accession.startswith(self.acc_prefix) and len(accession) == 11:
                        #file_name = dxencode.file_path_from_fid(fid,projectToo=True).split(':')[-1]
                        #print " - Already posted: " + file_name + " to " + accession
                        continue
            # No accession so try to match in encoded by submit_file_name and size
            #f_obj =  self.find_in_encode(fid,verbose)
            f_obj =  self.enc_file_find_find_by_dxid(fid)
            if f_obj == None:
                needed.append( (out_type,rep_tech,fid) )
            else:
                self.found[fid] = f_obj

        if verbose:
            print "Needed files:"
            print json.dumps(needed,indent=4)
        return needed

    def find_sw_versions(self,dxFile,dx_app=False,verbose=False):
        '''
        Finds the software versions associated with a file.
        Returns  { "software_versions": [ { "software": "star", "version": "2.4.0k" }, ... ] }
        '''
        sw_versions = {}
        # looks first in dx file property.
        # "SW" = {"DX applet": {"align-bwa-se.sh": "0.1.0"}, "samtools": "0.2.0", "bwa": "0.7.7-r441"}
        SW =  dxencode.dx_file_get_property("SW",None,dxfile=dxFile,return_json=True)
        if SW != None:
            if dx_app:
                if "DX applet" in SW:
                    SW = SW["DX applet"]
                else:
                    SW = {}
            versions = []
            for key in SW.keys():
                versions.append({"software": key, "version": SW[key]})
            sw_versions = { "software_versions": versions }
        else:
            # if no 'SW' property then try grepping from the log.
            regoop = '\* (\S+)\s+version:\s+(\S+)'
            if dx_app:
                regoop = '\* Running:\s(\S+)\s+\S*\[(\S+)\]'
            sw_versions = dxencode.get_sw_from_log(dxFile, regoop) # * STAR version: 2.4.0k
            
        if verbose:
            print "sw_versions: "
            #print dxFile
            print json.dumps(sw_versions,indent=4)
        return sw_versions

    def find_app_version(self,dxFile,verbose=False):
        '''
        Finds the app versions associated with a file.
        Returns  { "software": "align-bwa-se.sh", "version": "1.0.2"}
        '''
        app_version = {}
        sw_versions = self.find_sw_versions(dxFile,dx_app=True,verbose=verbose)
        if sw_versions != None:
            app_version =  sw_versions["software_versions"][0]
        if verbose:
            print "app_version: "
            print json.dumps(app_version,indent=4)
        return app_version
        
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
                acc_key = dxencode.dx_property_accesion_key(self.server)
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
        obj['lab'] = '/labs/encode-processing-pipeline/' # self.exp['lab']['@id'] Now hard-coded
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


    def get_qc_metrics(self,fid,job,verbose=False):
        '''Returns an object containing 'QC_metrics' info found in the file and or job.'''
        qc_metrics = None
        qc_for_job = None
        
        # if file level qc_metrics exist, then use those
        qc_for_file = dxencode.dx_file_get_property("QC",fid,return_json=True)
        if not qc_for_file:
            if "metadata" in job["output"]:
                try:
                    qc_for_job = json.loads(job["output"]["metadata"]) # NOTE: output class=string name='metadata' in json format
                except:
                    try:
                        qc_for_job = json.loads("{"+job["output"]["metadata"]+"}")
                    except:
                        print job["output"]["metadata"]
                        sys.exit(1)
            # if metadata is SW versions, then throw it out
            if qc_for_job != None and "DX applet" in qc_for_job:
                qc_for_job = None
                
        # We used to combine, giving qc_for_file priority.  Now we just return for qc_for_file if it exists.
        #if qc_for_job:
        #    qc_metrics = qc_for_job
        #    if qc_for_file:
        #        qc_metrics.update(qc_for_file)  # Note that qc_for_file take priority over qc_for_job
        #else:
        #    qc_metrics = qc_for_file
        if qc_for_file:
            qc_metrics = qc_for_file
        else:
            qc_metrics = qc_for_job
                    
        if qc_metrics and verbose:
            print "qc_metrics:"
            print json.dumps(qc_metrics,indent=4)
            
        return qc_metrics
        

    def pipeline_qualifiers(self,rep_tech,app_name=None):
        '''Determines and pipeline_qualifiers in ugly special-case code.'''
        if "exp" in self.obj_cache and rep_tech in self.obj_cache["exp"] \
        and "pipe_qualifiers" in self.obj_cache["exp"][rep_tech]:
            return self.obj_cache["exp"][rep_tech]["pipe_qualifiers"]
        #if "exp" in self.obj_cache and pipe_qualifiers in self.obj_cache["exp"]:
        #    return self.obj_cache["exp"]["pipe_qualifiers"]
            
        pipe_qualifiers = {}
        pipe_qualifiers["version"] = "1"
        if self.exp_type.startswith('long-rna-seq'): # FIXME: ugly special case
            pipe_qualifiers["version"] = "2"

        pipe_qualifiers["qualifier"] = ''
        if self.exp_type.startswith('long-rna-seq'): # FIXME: ugly special case
            #pipe_qualifiers["qualifier"] = '-unknown'
            if app_name and len(app_name):
                if app_name.endswith('-pe'):
                    pipe_qualifiers["qualifier"] = '-pe'
                elif app_name.endswith('-se'):
                    pipe_qualifiers["qualifier"] = '-se'
            if len(pipe_qualifiers["qualifier"]) == 0 and rep_tech and len(rep_tech):
                if rep_tech.startswith("rep"):
                    br_tr = rep_tech[3:]
                    (br,tr) = br_tr.split('_')
                    full_mapping = dxencode.get_full_mapping(self.exp_id,self.exp)
                    mapping = dxencode.get_replicate_mapping(self.exp_id,int(br),int(tr),full_mapping)
                    if "paired_ended" in mapping:
                        pipe_qualifiers["qualifier"] = '-pe'
                    else:
                        pipe_qualifiers["qualifier"] = '-se'
                else: #if rep_tech == "combined":
                    # FIXME: What about combined replicate only files???  No problem SO FAR since lrna has no combined steps
                    # Desperate attempt, but if the replicate level files are uploaded first then could try to look them up
                    if "exp" in self.obj_cache:
                        for try_rep in  ["rep1_1","rep2_1","rep1_2","rep2_2"]:
                            if try_rep in self.obj_cache["exp"] and "pipe_qualifiers" in self.obj_cache["exp"][try_rep]:
                                pipe_qualifiers["qualifier"] = self.obj_cache["exp"][try_rep]["pipe_qualifiers"]["qualifier"]
            if len(pipe_qualifiers["qualifier"]) == 0:
                print "Error: Can't determine pipeline qualifiers for '%s' replicate of '%' experiment" % \
                                                                                    (replicate,self.exp_type)
                sys.exit(1) 
                            
        if "exp" not in self.obj_cache:
             self.obj_cache["exp"] = {}
        if rep_tech not in self.obj_cache["exp"]:
            self.obj_cache["exp"][rep_tech] = {}
        self.obj_cache["exp"][rep_tech]["pipe_qualifiers"] = pipe_qualifiers
        #self.obj_cache["exp"]["pipe_qualifiers"] = pipe_qualifiers
        return pipe_qualifiers


    def enc_pipeline_find(self,dx_pipe_name,rep_tech):
        '''Finds a 'pipeline' encoded object.'''
        pipe = None
        # First lookup pipe_qualifiers:
        pipe_qualifiers = self.pipeline_qualifiers(rep_tech) 
        
        # Lookup pipeline as: /pipelines/encode:long-rna-seq-se-pipeline-2
        pipe_name = dx_pipe_name + pipe_qualifiers["qualifier"] + '-pipeline-' + pipe_qualifiers["version"].replace('.','-')
        pipe_alias = 'encode:' + pipe_name
        if pipe_alias in self.obj_cache:
            return self.obj_cache[pipe_alias]
        pipe = self.enc_lookup_json( 'pipelines/' + pipe_alias,must_find=True)
        if pipe:
            self.obj_cache[pipe_alias] = pipe
            print "  - Found pipeline: '%s'" % pipe_alias
        return pipe
       

    def enc_analysis_step_find(self,dx_app_name,dx_app_ver,dx_app_id,pipe_name):
        '''Finds the 'analysis_step' encoded object used in creating the file.'''
        ana_step = None
        ana_step_name = dx_app_name + '-v-' + dx_app_ver.replace('.','-')
        # NOTE: the special dnanexus: alias.  Because dx_ids will differ between projects, and the identical app can be
        #       rebuilt, the dnanexus: alias will instead be the app_name and app_version!
        #       The ana_step_name, on the otherhand, may be sanitized ('prep_star' => 'index_star')
        # /analysis-steps/unstranded-signal-star-v-1-0-1/ or /analysis-steps/dnanexus:bam-to-bigwig-unstranded-v-1-0-1/
        ana_step_alias = 'dnanexus:' + ana_step_name
        if ana_step_alias in self.obj_cache:
            return self.obj_cache[ana_step_alias]
        ana_step = self.enc_lookup_json( 'analysis-steps/' + ana_step_name,must_find=False)
        if not ana_step:
            ana_step = self.enc_lookup_json( 'analysis-steps/' + ana_step_alias,must_find=True)
        if ana_step:
            self.obj_cache[ana_step_alias] = ana_step
            print "  - Found ananlsys_step: '%s'" % ana_step_alias
        # Could attempt to create analysis-step
        #if not ana_step:  
        #    ana_step = {}
        #    ana_step['aliases'] = [ ana_step_alias  ]
        #    ana_step['name'] = ana_step_name  
        #    ana_step['pipeline'] = "/pipelines/encode:" + pipe_name
        #    #ana_step["pipeline"] = pipe
        #
        #    # Shoe-horn in:
        #    notes = {}
        #    notes['dx_app_name'] = dx_app_name
        #    notes['dx_app_ver'] = dx_app_ver
        #    ana_step["notes"] = json.dumps(notes)
        #    #    ana_step = dxencode.encoded_post_obj('analysis-step-run',ana_step, self.server, self.authid, self.authpw)
        return ana_step


    def enc_workflow_run_find_or_create(self,dx_wfr_id,rep_tech,test=False,verbose=False):
        '''Finds or creates the 'workflow_run' encoded object that actually created the file.'''
        wf_run = None
        wf_alias = 'dnanexus:' + dx_wfr_id
        # Find if you can:  /workflow-runs/dnanexus:analysis-BZKpZBQ0F1GJ7Z4ky8vBZpBx
        if "exp" in self.obj_cache and wf_alias in self.obj_cache["exp"]:
            wf_run = self.obj_cache["exp"][wf_alias]
        else:
            wf_run = self.enc_lookup_json( '/workflow-runs/' + wf_alias,must_find=False)
            if wf_run:
                self.obj_cache["exp"][wf_alias] = wf_run
                print "  - Found workflow_run: '%s'" % wf_alias
        
        if wf_run == None:
            wf_run = {}
            
            # Lookup pipeline as: /pipelines/encode:long-rna-seq-se-pipeline-2
            pipe = self.enc_pipeline_find(self.exp_type,rep_tech)
            
            wf_run['aliases'] = [ wf_alias ]
            wf_run['status'] = "finished"
            wf_run['pipeline'] = "/pipelines/encode:" + pipe["name"]
            wf_run["dx_analysis_id"] = dx_wfr_id
            wf_run["dx_workflow_id"] = dx_wfr_id # TODO: remove this when schema is updated
            # wf_run["software_version"] NO
            # wf_run["input_files"] NO
            # wf_run["dx_analysis_id"] NO
            # wf_run["dx_workflow_json"] NO

            # Look up the actual analysis
            dx_wf = dxpy.api.analysis_describe(dx_wfr_id)
            if dx_wf:  # Note that wf is created immediately before running, then last modified by dx to set status 'done'
                if "created" in dx_wf:
                    then = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(dx_wf.get("created")/1000.0))
                    wf_run["started_running"] = then
                if "modified" in dx_wf:
                    then = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(dx_wf.get("modified")/1000.0))
                    wf_run["stopped_running"] = then
            
            # shoe-horn into notes:
            notes = {}
            notes["notes_version"] = "1"
            notes['pipeline_name'] = pipe["name"]
            notes["dx_project_id"] = self.proj_id
            notes["dx_project_name"] = self.proj_name
            if dx_wf:
                notes["dx_cost"] = "$" + str(round(dx_wf.get("totalPrice"),2)) 
            wf_run["notes"] = json.dumps(notes)
            
            # Now post this new object
            self.obj_cache["exp"][wf_alias] = wf_run
            if test:
                print "  - Would post workflow_run: '%s'" % wf_alias
            else:
                try:
                    #wf_run["@type"] = ["item", "workflow_run"]
                    posted_wf_run = dxencode.encoded_post_obj('workflow_run',wf_run, self.server, self.authid, self.authpw)
                except:
                    print "Failed to post workflow_run: '%s'" % wf_alias
                    sys.exit(1)
                print "  - Posted workflow_run: '%s'" % wf_alias
                #sys.exit(1)  # TEST one case!
            
        if verbose:
            print "ENCODEd 'workflow_run':"
            print json.dumps(wf_run,indent=4)
            if "notes" in wf_run:
                wf_run_notes = json.loads(wf_run.get("notes"))
                print "ENCODEd 'workflow_run[notes]':"
                print json.dumps(wf_run_notes,indent=4)
        return wf_run
        
        
    def enc_step_run_find_or_create(self,job,dxFile,rep_tech,test=False,verbose=False):
        '''Finds or creates the 'analysis_step_run' encoded object that actually created the file.'''
        step_run = None
        job_id = job.get('id')
        step_alias = 'dnanexus:' + job_id
        if "exp" in self.obj_cache and step_alias in self.obj_cache["exp"]:
            step_run = self.obj_cache["exp"][step_alias]
        else:
            step_run = self.enc_lookup_json( '/analysis-step-runs/' + step_alias,must_find=False)
            if step_run:
                if "exp" not in self.obj_cache:
                     self.obj_cache["exp"] = {}
                self.obj_cache["exp"][step_alias] = step_run
                print "  - Found step_run: '%s'" % step_alias

        if step_run == None:  
            step_run = {}
            step_run['aliases'] = [ step_alias ]
            step_run['status'] = "finished"
            dx_app_name = job.get('executableName')
            
            # MUST determine special case pipeline qualifiers before proceeding!
            pipe_qualifiers = self.pipeline_qualifiers(rep_tech,dx_app_name) 
            # Note that qualifiers are not used here, but will be in workflow_run creation

            # Find or create the workflow
            dx_wfr_id = job.get('analysis')
            #wf_run = self.enc_workflow_run_find_or_create(dx_wfr_id,rep_tech,test=False,verbose=verbose)
            wf_run = self.enc_workflow_run_find_or_create(dx_wfr_id,rep_tech,test=self.test,verbose=verbose)
            wf_run_notes = json.loads(wf_run["notes"])
            step_run["workflow_run"] = "/workflow-runs/dnanexus:" + dx_wfr_id
            
            # Find analysis_step
            dx_app_id = job.get('applet')
            #dx_app_ver = dxencode.get_sw_from_log(dxFile, '\* Running:\s(\S+)\s+\S*\[(\S+)\]') # * Running: align-star-pe.sh [v2.0.0]
            dx_app_ver = self.find_app_version(dxFile)
            if dx_app_ver and 'version' in dx_app_ver:
                dx_app_ver = str( dx_app_ver.get('version') )
                if dx_app_ver[0] == 'v':
                    dx_app_ver = dx_app_ver[1:]
            if not dx_app_ver or not isinstance(dx_app_ver, str) or len(dx_app_ver) == 0:
                print "ERROR: cannot find applet version %s in the log" % ( type(dx_app_ver) )
                sys.exit(0)
            ana_step = self.enc_analysis_step_find(dx_app_name,dx_app_ver,dx_app_id,wf_run_notes["pipeline_name"])
            step_run['analysis_step'] = "/analysis-steps/" + ana_step['name']
            
            # applet details:
            applet_details = {}
            #step_run["dx_applet_details"] = {}
            applet_details["dx_job_id"] = job_id
            then = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(job.get('startedRunning')/1000.0))
            applet_details["started_running"] = then
            then = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(job.get('stoppedRunning')/1000.0))
            applet_details["stopped_running"] = then
            applet_details["dx_status"] = "finished"

            # Parameters:
            params = {}
            inputs = job.get("originalInput")
            for name in inputs.keys():
                #if isinstance(inputs[name],str) or isinstance(inputs[name],int) or isinstance(inputs[name],unicode):
                if not isinstance(inputs[name],dict):
                    params[name] = inputs[name]
            if len(params) > 0:
                applet_details["parameters"] = params
                
            # applet json?
            #dx_app = dxpy.api.applet_describe(dx_app_id)
            #if dx_app:  # Note that wf is created immediately before runiing, then last modified by dx to set status
            #    applet_details["dx_app_json"] = dx_app
            step_run["dx_applet_details"] = [ applet_details ]
            
            # shoe-horn into notes:
            notes = {}
            notes["notes_version"] = "1"
            notes["dx_app_name"] = dx_app_name
            notes['dx_app_version'] = dx_app_ver
            notes["dx_app_id"] = dx_app_id
            notes["step_name"] = ana_step['name']
            notes['pipeline_name'] = wf_run_notes["pipeline_name"]
            notes["dx_analysis_id"] = dx_wfr_id
            notes["dx_project_id"] = self.proj_id
            notes["dx_project_name"] = self.proj_name
            notes["dx_cost"] = "$" + str(round(job['totalPrice'],2))
            step_run["notes"] = json.dumps(notes)
            
            # Now post this new object
            self.obj_cache["exp"][step_alias] = step_run
            if test:
                print "  - Would post step_run: '%s'" % step_alias
            else:
                try:
                    #step_run["@type"] = ["item", "analysis_step_run"]
                    posted_step_run = dxencode.encoded_post_obj('analysis_step_run',step_run, self.server, self.authid, self.authpw)
                except:
                    print "Failed to post step_run: '%s'" % step_alias
                    sys.exit(1)
                print "  - Posted step_run: '%s'" % step_alias
                #sys.exit(1)  # TEST one case!

        if step_run and verbose:
            print "ENCODEd 'analysis_step_run':"
            print json.dumps(step_run,indent=4)
            if "notes" in step_run:
                step_run_notes = json.loads(step_run.get("notes"))
                print "ENCODEd 'analysis_step_run[notes]':"
                print json.dumps(step_run_notes,indent=4)
        return step_run
        
        
    def make_payload_obj(self,out_type,rep_tech,fid,verbose=False):
        '''Returns an object for submitting a file to encode, with all dx info filled in.'''
        payload = {}
        payload['dataset'] = self.exp_id
        payload["output_type"] = out_type

        dx_obj = dxencode.description_from_fid(fid)

        # get job for this file
        job = dxencode.job_from_fid(fid)

        payload["file_format"] = self.file_format(dx_obj["name"])
        if payload["file_format"] == None:
            print "Warning: file %s has unknown file format!" % dxencode.file_path_from_fid(fid)
        payload["derived_from"] = self.find_derived_from(fid,job, verbose)
        payload['submitted_file_name'] = dxencode.file_path_from_fid(fid,projectToo=True)
        payload['file_size'] = dx_obj["size"]
        #payload['md5sum'] = calculated_md5 # Done i  validate_post applet
        if self.genome == None:
            print "Warning: could not determine genome assembly! Add properties to reference files."
            #sys.exit(1)
        else:
            payload['assembly'] = self.genome
        if self.annotation != None:
            payload['genome_annotation'] = self.annotation

        dxFile = dxencode.file_handler_from_fid(fid)
        #versions = dxencode.get_sw_from_log(dxFile, '\* (\S+)\s+version:\s+(\S+)') # * STAR version: 2.4.0k
        versions = self.find_sw_versions(dxFile)
        notes = dxencode.create_notes(dxFile, versions)
        notes["notes_version"] = "5" # Cricket requests starting at "5", since earlier files uploads were distingusihed by user
        if 'totalPrice' in job:
            notes['dx_cost'] = "$" + str(round(job['totalPrice'],2))

        # Put QC metrics in notes
        qc_metrics = self.get_qc_metrics(fid,job)
        if qc_metrics:
            notes['QC_metrics'] = qc_metrics
        
        #print "  - Adding encoded information."
        payload = self.add_encoded_info(payload,rep_tech,fid)

        # Find or create step_run object
        step_run = self.enc_step_run_find_or_create(job,dxFile,rep_tech,test=self.test,verbose=verbose)
        #step_run = self.enc_step_run_find_or_create(job,dxFile,rep_tech,test=False,verbose=verbose)
        if step_run:
            if "dx_applet_details" in step_run:
                payload["step_run"] = "/analysis-step-runs/dnanexus:" + step_run["dx_applet_details"][0].get("dx_job_id")
            if "notes" in step_run:
                step_run_notes = json.loads(step_run.get("notes"))
                # analysis_step and pipeline are calculated properties.
                #payload["analysis_step"] = "/analysis-steps/" + step_run_notes.get("step_name")
                #payload["pipeline"] = "/pipelines/encode:" + step_run_notes.get("pipeline_name")
                if "dx_analysis_id" in step_run_notes:
                    notes["workflow_run"] = "/workflow-runs/dnanexus:" + step_run_notes.get("dx_analysis_id")
                elif "dx_workflow_id" in step_run_notes:
                    notes["workflow_run"] = "/workflow-runs/dnanexus:" + step_run_notes.get("dx_workflow_id")

        notes["dx_project_id"] = self.proj_id
        notes["dx_project_name"] = self.proj_name
        payload['notes'] = json.dumps(notes)
        
        if verbose:
            print "payload:"
            print json.dumps(payload,indent=4)
            print "payload[notes]:"
            print json.dumps(notes,indent=4)
        return payload


    def file_mark_accession(self,fid,accession,test=True):
        '''Adds/replaces accession to a file's properties.'''
        acc_key = dxencode.dx_property_accesion_key(self.server)
        path = dxencode.file_path_from_fid(fid)
        acc = dxencode.dx_file_set_property(fid,acc_key,accession,add_only=True,test=test)
        if acc == None or acc != accession:
            print "Error: failed to update %s for file %s to '%s'" % (path,acc_key,accession)
        elif test:
            print "  - Test flag %s with %s='%s'" % (path,acc_key,accession)
        else:
            print "  - Flagged   %s with %s='%s'" % (path,acc_key,accession)


    def can_skip_validation(self,exp_type,output_type,test=True):
        '''Returns True if validation can be skipped.'''
        if exp_type.startswith("long-rna-seq"):
            return True
        return (output_type in self.SKIP_VALIDATE)

    def file_post(self,fid,payload,test=True):
        '''Posts a file to encoded.'''
        path = payload['submitted_file_name'].split(':')[1]
        derived_count = len(payload["derived_from"])
        skip_validate = self.can_skip_validation(self.exp_type,payload["output_type"])
        job_name = "Post "+path+" to "+self.server_key
        if not skip_validate:
            job_name = job_name + " (must validate)"
        if test:
            print "  - Test %s (derived:%d)" % (job_name,derived_count)
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
                folder=out_folder,name=job_name)
            print "  - Job: %s <%s> (derived:%d)" % \
                                                     (job_name,job.id, derived_count)
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
        self.ignore = False
        if args.ignore_properties:
            print "Ignoring DXFile properties (will post to test server)"
            self.ignore = args.ignore_properties
            self.server_key = 'test' # mandated because option is dangerous
            
        self.server_key = args.server
        self.authid, self.authpw, self.server = dxencode.processkey(self.server_key)
        
        if self.server_key == "www":
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
            self.exp_id = exp_id
            self.obj_cache["exp"] = {}  # clear exp cache, which will hold exp specific wf_run and step_run objects
            # 1) Lookup experiment type from encoded, based on accession
            print "Working on %s..." % self.exp_id
            self.exp = dxencode.get_exp(self.exp_id,must_find=False,key=self.server_key)
            if self.exp == None or self.exp["status"] == "error":
                print "Unable to locate experiment %s in encoded" % self.exp_id
                continue
            self.exp_type = dxencode.get_exp_type(self.exp_id,self.exp,self.EXPERIMENT_TYPES_SUPPORTED)
            if self.exp_type == None:
                continue

            # 2) Locate the experiment accession named folder
            # NOTE: genome and annotation are not known for this exp yet, so the umbrella folder is just based on exp_type
            self.umbrella_folder = dxencode.umbrella_folder(args.folder,self.FOLDER_DEFAULT,self.proj_name,self.exp_type)
            self.exp_folder = dxencode.find_exp_folder(self.project,exp_id,self.umbrella_folder,warn=True)
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
            files_expected = self.find_expected_files(self.exp_folder, self.replicates, verbose=args.verbose)
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
                elif accession == "NOT POSTED":
                    print "* HALTING %s - validation failure prevented posting and could compromise 'derived_from'" % \
                                                                                    (self.exp_id)
                    halted += 1
                    break

                # d) Update dnanexus file with file accession tag.
                if not args.test:
                    post_count += 1
                    
                self.file_mark_accession(fid,accession,args.test)  # This should have already been set by validate_post

                #if file_count >= 1:  # Short circuit for test
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

