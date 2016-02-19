#!/usr/bin/env python2.7
# splashdown.py 1.0.0
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
from base64 import b64encode
import commands

import dxpy
import dxencode

class Splashdown(object):
    '''
    Splashdown module posts from dnanexus to ENCODEd,  all files available and necessry for
    a given experiment .
    '''
    TOOL_IS = 'splashdown'
    HELP_BANNER = "Handles splashdown of launched pipeline runs for supported experiment types. " + \
                  "Can be run repeatedly and will only try to post result files that have not been previously posted. " 
    ''' This help banner is displayed by get_args.'''

    SERVER_DEFAULT = 'test'
    '''This the default server to post files to.'''

    FOLDER_DEFAULT = "/"
    '''Where to start the search for experiment folders.'''

    EXPERIMENT_TYPES_SUPPORTED = [ 'long-rna-seq', 'small-rna-seq', 'rampage', 'dna-me' ] #, 'dnase' ]
    '''This module supports only these experiment (pipeline) types.'''

    SKIP_VALIDATE = {"transcription start sites":'bed'}
    '''Some output_types cannot currently be validated, theoretically'''

    # Pipeline specifications include order of steps, steps per replicate, combined steps and
    # within steps, the output_type: file_glob that define expected results.
    # Note: that some steps have multiple files with the same output_type (e.g. hotspot: bed & bb). 
    #       When this happens, key on "output_type|format|format_type": file_glob 
    #       (e.g. "hotspot|bed|narrowPeak": "*_hotspot.bed" and "hotspot|bb|narrowPeak": "*_hotspot.bb")
    #       TODO: This could be done more gracefully.
    PIPELINE_SPECS = {
         "long-rna-seq": {
            "step-order": [ "align-tophat","signals-top-se","signals-top-pe",
                            "align-star","signals-star-se","signals-star-pe","quant-rsem","mad-qc"],
            "replicate":  {
                "align-tophat":    { "alignments":                                "*_tophat.bam"                 },
                "signals-top-se":  { "signal of all reads":                       "*_tophat_all.bw",
                                     "signal of unique reads":                    "*_tophat_uniq.bw"             },
                "signals-top-pe":  { "minus strand signal of all reads":          "*_tophat_minusAll.bw",
                                     "plus strand signal of all reads":           "*_tophat_plusAll.bw",
                                     "minus strand signal of unique reads":       "*_tophat_minusUniq.bw",
                                     "plus strand signal of unique reads":        "*_tophat_plusUniq.bw"         },
                "signals-star-se": { "signal of all reads":                       "*_star_genome_all.bw",
                                     "signal of unique reads":                    "*_star_genome_uniq.bw"        },
                "signals-star-pe": { "minus strand signal of all reads":          "*_star_genome_minusAll.bw",
                                     "plus strand signal of all reads":           "*_star_genome_plusAll.bw",
                                     "minus strand signal of unique reads":       "*_star_genome_minusUniq.bw",
                                     "plus strand signal of unique reads":        "*_star_genome_plusUniq.bw"    },
                "align-star":      { "alignments":                                "*_star_genome.bam",
                                     "transcriptome alignments":                  "*_star_anno.bam"              },
                "quant-rsem":      { "gene quantifications":                      "*_rsem.genes.results",
                                     "transcript quantifications":                "*_rsem.isoforms.results"      }  },
            "combined":   {
                "mad-qc":          { "QC_only":                                   "*_mad_plot.png"               }  },
        },
        "small-rna-seq": {
            "step-order": [ "align","signals","mad_qc"],
            "replicate":  {
                "align":           { "alignments":                                "*_srna_star.bam",
                                     "gene quantifications":                      "*_srna_star_quant.tsv"        },
                "signals":         { "plus strand signal of all reads":           "*_small_plusAll.bw",
                                     "minus strand signal of all reads":          "*_small_minusAll.bw",
                                     "plus strand signal of unique reads":        "*_small_plusUniq.bw",
                                     "minus strand signal of unique reads":       "*_small_minusUniq.bw"         }  },
            "combined":   { 
                "mad_qc":          { "QC_only":                                   "*_mad_plot.png"               }  },
        },
        "rampage": {
            "step-order": [ "align","signals","peaks","idr","mad_qc"],
            "replicate":  {
                "align":           { "alignments":                                "*_rampage_star_marked.bam" },
                "signals":         { "plus strand signal of all reads":           "*_rampage_5p_plusAll.bw",
                                     "minus strand signal of all reads":          "*_rampage_5p_minusAll.bw",
                                     "plus strand signal of unique reads":        "*_rampage_5p_plusUniq.bw",
                                     "minus strand signal of unique reads":       "*_rampage_5p_minusUniq.bw" },
                "peaks":           { "transcription start sites|gff|gff3":        "*_rampage_peaks.gff",
                                     "transcription start sites|bed|tss_peak":    "*_rampage_peaks.bed",
                                     "transcription start sites|bigBed|tss_peak": "*_rampage_peaks.bb",
                                     "gene quantifications":                      "*_rampage_peaks_quant.tsv" } },
            "combined":   {
                "idr":             { "transcription start sites|bed|idr_peak":    "*_rampage_idr.bed",
                                     "transcription start sites|bigBed|idr_peak": "*_rampage_idr.bb" },
                "mad_qc":          { "QC_only":                                   "*_rampage_mad_plot.png" }  },
        },
        "dna-me": {
            "step-order": [ "align","quantification","corr"], # How to: 1) combine 3 steps into 1; 2) tech lvl, bio lvl, exp lvl
            "replicate":  {
                "align":           { "alignments":                                "*_bismark.bam" }, # *may* have samtools_flagstat, samtools_stats, Don't wan't bismark_map
                "quantification":  { "methylation state at CpG|bigBed|bedMethyl": "*_bismark_biorep_CpG.bb",      # All have: samtools_flagstat, bismark_map
                                     "methylation state at CpG|bed|bedMethyl":    "*_bismark_biorep_CpG.bed.gz",  # All have: samtools_flagstat, bismark_map
                                     "methylation state at CHG|bigBed|bedMethyl": "*_bismark_biorep_CHG.bb",      # All have: samtools_flagstat, bismark_map
                                     "methylation state at CHG|bed|bedMethyl":    "*_bismark_biorep_CHG.bed.gz",  # All have: samtools_flagstat, bismark_map
                                     "methylation state at CHH|bigBed|bedMethyl": "*_bismark_biorep_CHH.bb",      # All have: samtools_flagstat, bismark_map
                                     "methylation state at CHH|bed|bedMethyl":    "*_bismark_biorep_CHH.bed.gz",  # All have: samtools_flagstat, bismark_map
                                     "signal":                                    "*_bismark_biorep.bw" } },      # All have: samtools_flagstat, bismark_map
            "combined":   {
                "corr":            { "QC_only":                                   "*_CpG_corr.txt" }  }, # Not yet defined in encodeD
        },
    }
    
    # Step children are steps that should be combined with their parent step rather than be treated as a separate job 
    STEP_CHILDREN = {
        "dme-cx-to-bed":        "dme-extract-pe",
        "dme-cx-to-bed-alt":    "dme-extract-se",
        "dme-bg-to-signal":     "dme-extract-pe",
        "dme-bg-to-signal-alt": "dme-extract-se",
    } 

    ASSEMBLIES_SUPPORTED = { "hg19": "hg19", "GRCh38": "GRCh38", "mm10": "mm10" }
    '''This module supports only these assemblies.'''

    ANNOTATIONS_SUPPORTED = [ 'V24', 'V19', 'M2', 'M3', 'M4' ]
    '''This module supports only these annotations.'''
    
    REQUIRE_ANNOTATION = [ 'long-rna-seq','small-rna-seq','rampage' ]
    '''These assays require an annotation.'''

    FORMATS_SUPPORTED = ["bam","bed","bigBed","bigWig","fasta","fastq","gff","gtf","hdf5","idat","rcc","CEL",
                         "tsv","csv","sam","tar","wig"]
    FORMAT_TYPES_SUPPORTED = ["bed6","bed9","bed12","bedGraph","bedLogR","bedMethyl","broadPeak",
                        "enhancerAssay","gappedPeak","gff2","gff3","narrowPeak" ]
    # TODO: format types need to support rampage: tss, and idr
    EXTENSION_TO_FORMAT = { 
        "2bit":       "2bit",
        "cel.gz":     "CEL",
        "bam":        "bam",
        "bed.gz":     "bed",     "bed":     "bed",
        "bigBed":     "bigBed",  "bb":      "bigBed",
        "bigWig":     "bigWig",  "bw":      "bigWig",
        "csfasta.gz": "csfasta",
        "csqual.gz":  "csqual",
        "fasta.gz":   "fasta",   "fa.gz":   "fasta",  "fa": "fasta",
        "fastq.gz":   "fastq",   "fq.gz":   "fastq",  "fq": "fastq",
        "gff.gz":     "gff",     "gff":     "gff",
        "gtf.gz":     "gtf",     "gtf":     "gtf",
        "h5":         "hdf5",
        "idat":       "idat",
        "rcc":        "rcc",
        "tar.gz":     "tar",     "tgz":     "tar",
        "tsv":        "tsv",     "results": "tsv",
        "csv":        "csv",
        "wig.gz":     "wig",     "wig":     "wig",
        "sam.gz":     "sam",     "sam":     "sam"
        }
    #EXTENSION_TO_FORMAT = { "bb":"bigBed", "bw":"bigWig",
    #                        "fa":"fasta","fq":"fastq","results":"tsv",
    #                        "gff": "gtf" }
    '''List of supported formats, and means of recognizing with file extensions.'''
    
    # Each QC object has its own particulars:
    # The key should match a key found in a file details json
    # type: (optional) the collection of quality_metric objects in encodeD. Collection will be type+'_quality_metric'. 
    #       Default: lower(key) 
    # files: (optional) object containing "inputs" and/or "results". 
    #        "results" is the usual case where the qc_obj is in the dx file 'details'. 
    #        "inputs" points to a list of named inputs to the job. 
    # blob: (optional) The 'pattern' to match and (optional) 'mime_ext' of a file to determine mime_type.
    # singleton: (optional) There is only one value in dx detail json.  It will be named as such.
    # include: (optional) A list of the only props in dx json to be placed in enc object. Default: include all
    # exclude: (optional) A list of props in dx json but not to be placed in enc object.  Default: exclude none
    # props: (optional) Select and redefine the names of any properties in dx details object that should be in enc object.
    #        Default: include all from dx details json and use the same property names
    QC_SUPPORTED = {
        "STAR_log_final":       { 
                                    "type":"star", 
                                    "files": { "results": "detail"}, 
                                    "blob":  { "pattern": "/*_star_Log.final.out",  
                                               "mime_ext": "txt" }, 
                                    "exclude": [ "Finished on","Started job on","Started mapping on"],
                                },
        "MAD.R":                { 
                                    "type":"mad",
                                    "files": { "inputs": ["quants_a", "quants_b"] },        
                                    "blob": { "pattern": "/*_mad_plot.png" }, 
                                },
        "IDR_summary":          { 
                                    "files": {"results": "detail", "inputs": ["peaks_a", "peaks_b"] },        
                                    "blob": { "pattern": "/*_idr.png" }, 
                                },
        # TODO: characterized DNase qc_metrics better.  Especially files.
        "samtools_flagstats":   { 
                                    "files": {"results": "detail"}, 
                                    "blob": { "pattern": "/*_qc.txt" } 
                                },
        "tophat_flagstat":      { 
                                    "type":"samtools_flagstats",
                                    "only_for": "_tophat.bam",
                                    "files": {"results": "detail"}, 
                                    "blob": { "pattern": "/*_tophat_flagstat.txt" }, 
                                },
        "star_genome_flagstat": { 
                                    "type":"samtools_flagstats",
                                    "only_for": "_star_genome.bam",
                                    "files": {"results": "detail"}, 
                                    "blob": { "pattern": "/*_star_genome_flagstat.txt" }, 
                                },
        "star_anno_flagstat":   { 
                                    "type":"samtools_flagstats",
                                    "only_for": "_star_anno.bam",
                                    "files": {"results": "detail"}, 
                                    "blob": { "pattern": "/*_star_anno_flagstat.txt" }, 
                                },
        "star_srna_flagstat": { 
                                    "type":"samtools_flagstats",
                                    "only_for": [ "_srna_star.bam", "_srna_star_quant.tsv" ],
                                    "files": {"results": "detail"}, 
                                    "blob": { "pattern": "/*_srna_star_flagstat.txt" }, 
                                },
        "star_rampage_flagstat": { 
                                    "type":"samtools_flagstats",
                                    "only_for": "_rampage_star_marked.bam",
                                    "files": {"results": "detail"}, 
                                    "blob": { "pattern": "/*_rampage_star_marked_flagstat.txt" }, 
                                },
        "samtools_stats":       { "files": {"results": "detail"}, "blob": { "pattern": "/*_qc.txt"          } },
        "bismark_map":          { 
                                    "type":"bismark",
                                    "files": {"results": "detail"}, 
                                    "blob": { "pattern": "/*_bismark_biorep_map_report.txt" },
                                    "include": [
                                        "C methylated in CHG context","lambda C methylated in CHG context", 
                                        "C methylated in CHH context","lambda C methylated in CHH context",
                                        "C methylated in CpG context","lambda C methylated in CpG context",
                                        "Mapping efficiency",         "lambda Mapping efficiency",
                                        "Sequences analysed in total","lambda Sequences analysed in total" ],
                                },
        "bedmethyl_corr":       { 
                                    "type":"cpg_correlation",
                                    "files": {"inputs": ["CpG_A", "CpG_B"]}, 
                                    "blob": { "pattern": "/*_CpG_corr.txt"},
                                },
        "hotspot":              { "files": {"results": "detail"}, "blob": { "pattern": "/*_hotspot_qc.txt"  } },
        "edwBamStats":          { "files": {"results": "detail"}, "blob": { "pattern": "/*_qc.txt"          } },
        "edwComparePeaks":      { "files": {"results": "detail"}, "blob": { "pattern": "/*_qc.txt"          },"exclude": [ "aFileName","bFileName"] },
        "bigWigCorrelate":      { "files": {"results": "detail"}, "blob": { "pattern": "/*_qc.txt"          }, "singleton":"Pearson's correlation" },
        "pbc":                  { "files": {"results": "detail"}, "blob": { "pattern": "/*_qc.txt"          } },
        "peak_counts":          { 
                                    "type":"dnase_peaks", 
                                    "files": {"results": "detail"},        
                                    "blob": { "pattern": "/*_qc.txt"          } 
                                },
        "phantompeaktools_spp": { 
                                    "files": {"results": "detail"}, 
                                    "blob": { "pattern": "/*_qc.txt" }, 
                                },
        # How to deal with *_qc.txt? file name stripped of ending add "_qc.txt" and see if it is found.
    }
    
    # blob attachments have mime_types recognizable by file extension
    EXTENSION_TO_MIME = { 
        "png":       "image/png",
        "txt":       "text/plain",
        "tsv":       "text/tab-separated-values",
        "pdf":       "application/pdf",
        "jpg":       "image/jpeg",
        "gif":       "image/gif",
        "tiff":      "image/tiff",
        "html":      "text/html",
        "as":        "text/autosql",
    }

    PRIMARY_INPUT_EXTENSION = [ "fastq","fq"]
    '''List of file extensions used to recognize primary inputs to parse accessions.'''
    
    REFERENCE_EXTENSIONS = [ "tgz'", "tar","sizes", "fa", "fasta", "gtf" ]
    '''List of file extensions used to recognize reference files  Could be a problem for gtf's!'''
    
    REFERERNCE_PROJECT_ID = 'project-BKbfvF00Qy5PZF5V7Gv003v7'
    '''Simple way to recognize if a reference file is in the right project, and therefore has accessions.'''
    
    APPEND_FLAG = "APPEND_ME"
    '''Special flag allows appending rather than replacing derived_by on recovery.'''
    
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
        self.workflow_runs_created = 0
        self.step_runs_created = 0
        self.way_back_machine = False # Don't support methods/expectations used on very old runs.  Only modern methods!      
        self.exp_files = None # Currently only used by 'recovery' and the way_back_machine
        self.found = {} # stores file objects from encode to avoid repeated lookups # TODO: place in obj_cache?
        logging.basicConfig(format='%(asctime)s  %(levelname)s: %(message)s')
        dxencode.logger = logging.getLogger(__name__ + '.dxe') # I need this to avoid some errors
        dxencode.logger.addHandler(logging.StreamHandler()) #logging.NullHandler)
        print



    def get_args(self,parse=True):
        '''Parse the input arguments.'''
        ### PIPELINE SPECIFIC
        ap = argparse.ArgumentParser(description=self.HELP_BANNER + "All results " +
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

        ap.add_argument('-g','--genome',
                        help="The genome assembly that files were aligned to (default: discovered if possible)",
                        default=None,
                        required=False)

        ap.add_argument('--test',
                        help='Test run only, do not launch anything.',
                        action='store_true',
                        required=False)

        ap.add_argument('--start_at',
                        help="Start processing with this file name (or possibly accession).",
                        default=None,
                        required=False)

        ap.add_argument('--files',
                        help="Just upload this number of files (default: all)",
                        type=int,
                        default=0,
                        required=False)

        ap.add_argument('-w','--way_back_machine',
                        help="Use the 'way back machine' to find files posted long ago.",
                        action='store_true',
                        required=False)

        ap.add_argument('--qc_from_file',
                        help="If available, generate qc from file instead of dx details.",
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

        if parse:
            return ap.parse_args()
        else:
            return ap


    def find_replicate_folders(self,exp_folder,verbose=False):
        '''Returns a sorted list of replicate folders which are sub-folders of an exp folder.'''
        # normalize
        try:
            sub_folders = self.project.list_folder(exp_folder)['folders']
        except:
            if verbose:
                print >> sys.stderr, "No subfolders found for %s" % exp_folder
            return []

        replicates = []
        for path in sorted(sub_folders):
            folder = path.split('/')[-1]
            if folder.startswith('rep'):
                if len(folder[3:].split('_')) == 2:
                    replicates.append( folder )
        if verbose:
            print >> sys.stderr, "Replicate folders:"
            print >> sys.stderr, json.dumps(replicates,indent=4)
        return replicates


    def pipeline_specification(self,args,exp_type,exp_folder,verbose=False):
        '''Sets the pipeline specification object for this experiment type.'''

        # Start with dict containing common variables
        #self.expected = copy.deepcopy(self.PIPELINE_SPECS[exp_type])

        pipeline_specs = self.PIPELINE_SPECS.get(exp_type)
        self.annotation = None  # TODO: if appropriate, need way to determine annotation

        if verbose:
            print >> sys.stderr, "Pipeline specification:"
            print >> sys.stderr, json.dumps(pipeline_specs,indent=4)
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


    def file_mime_type(self,file_name):
        '''Try to determine mime_type from file extension.'''
        ext = file_name.split(".")[-1]
        if ext in self.EXTENSION_TO_MIME.keys():
            return self.EXTENSION_TO_MIME[ext]
        return None


    def find_genome_annotation(self,posting_fid,derived_from_file_dict):
        '''Try to determine genome from input file properties.'''
        # TODO: currently done in derived_from which is only run on needed files
        #       much change to do this on expected files.
        properties = derived_from_file_dict["properties"]
        msg = ""
        if self.genome == None:
            if "genome" in properties:
                genome = properties["genome"]
                if genome in self.ASSEMBLIES_SUPPORTED.keys():
                    self.genome = self.ASSEMBLIES_SUPPORTED[genome]
            else:
                file_path_parts = dxencode.file_path_from_fid(posting_fid).split('/')
                for assembly in self.ASSEMBLIES_SUPPORTED.keys():
                    if assembly in file_path_parts:
                        self.genome = assembly
                        break
            if self.genome != None:
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
            if token != "QC_only": 
                if self.file_format(file_globs[token]) == None:
                    print "Error: file glob %s has unknown file format! Please fix" % file_globs[token]
                    sys.exit(1)
            if verbose:
                print >> sys.stderr, "-- Looking for %s" % (result_folder + file_globs[token])
            fid = dxencode.find_file(result_folder + file_globs[token],self.proj_id, recurse=False)
            if fid != None:
                QC_only = (token == "QC_only") # Use new qc_object posting methods  
                step_files.append( (token,rep_tech,fid,QC_only) )
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
            print >> sys.stderr, "Expected files:"
            print >> sys.stderr, json.dumps(expected,indent=4)
        return expected


    def enc_lookup_json(self, path,frame='object',must_find=False):
        '''Commonly used method to get a json object from encodeD.'''
        url = self.server + path + '/?format=json&frame=' + frame
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


    def enc_file_find_by_dxid(self,dx_fid):
        '''Finds a encoded 'file' object by dnanexus alias.'''
        file_obj = None
        
        file_alias = 'dnanexus:' + dx_fid
        file_obj = self.enc_lookup_json( 'files/' + file_alias,must_find=False)
        return file_obj
       

    def find_files_using_way_back_machine(self,fid,exp_files,verbose=False):
        '''Looks for file in encode with the same submitted file name as the fid has.'''
        file_path = dxencode.file_path_from_fid(fid,projectToo=True)
        file_size = dxencode.description_from_fid(fid).get('size')
        file_name = file_path.split('/')[-1]
        
        if verbose:
            print >> sys.stderr, "Looking for '%s' of size: %d" % (file_name, file_size)

        found_file = None
        for enc_file in exp_files:
            if enc_file.get("submitted_file_name").endswith(file_name):
                if enc_file.get('file_size') == file_size:
                    # Check file.notes['dx_id']
                    if 'notes' in enc_file:
                        notes = json.loads(enc_file["notes"])
                        if 'dx_id' not in notes or notes['dx_id'] == fid:
                            found_file = enc_file
                            break
                    else:
                        found_file = enc_file
                        break
            if verbose:
                print >> sys.stderr, "  %s %d" % (enc_file.get("submitted_file_name"),enc_file.get('file_size'))
                
        if found_file != None and verbose:
            print >> sys.stderr, "Found file:"
            print >> sys.stderr, json.dumps(found_file,indent=4)
        return found_file


    def enc_file_add_alias(self,fid,accession,f_obj,remove=False,test=True):
        '''Updates ENCODEd file with its 'fid' based alias.'''
        fid_alias = 'dnanexus:' + fid
        update_payload = {}

        if 'aliases' in f_obj:
            if fid_alias in f_obj['aliases'] and not remove:
                return False
            elif fid_alias not in f_obj['aliases'] and remove:
                return False
            update_payload['aliases'] = f_obj['aliases']
        else:
            update_payload['aliases'] = []
        if not remove:
            update_payload['aliases'].append(fid_alias)
        elif fid_alias in update_payload['aliases']:
            update_payload['aliases'].remove(fid_alias)
        if not test:
            ret = dxencode.encoded_patch_obj(accession, update_payload, self.server, self.authid, self.authpw)
            #if ret == accession:
            if not remove:
                print "  * Updated ENCODEd '"+accession+"' with alias "+fid_alias+"."
            else:
                print "  * Removed ENCODEd '"+accession+"' alias "+fid_alias+"."
        else:
            if not remove:
                print "  * Would update ENCODEd '"+accession+"' with alias "+fid_alias+"."
            else:
                print "  * Would remove ENCODEd '"+accession+"' alias "+fid_alias+"."
        return True


    def find_posted_files(self,files_expected,test=True,verbose=False):
        '''Returns the tuple list of files already posted to ENCODEd.'''
        #verbose=True
        # get all files associated with the experiment up front, just in case it is needed:
        self.exp_files = dxencode.get_enc_exp_files(self.exp,key=self.server_key)
        if len(self.exp_files) == 0:
            print "* ERROR: found no files associated with experiment: " + self.exp_id

        posted = []
        self.found = {}
        self.revoked = []
        for (out_type, rep_tech, fid, QC_only) in files_expected:
            if QC_only: # Use new qc_object posting methods 
                continue
            if fid in self.found.keys():
                continue  # No need to handle the same file twice
            if verbose:
                print >> sys.stderr, "* DEBUG Working on: " + dxencode.file_path_from_fid(fid,projectToo=True)
                
            # Posted files should have accession in properties
            # TODO: make use of new file_get_acession() insteead
            #accession = self.file_get_accession(qc_fid,verify=True,fake_if_needed=False)
            #if fid in self.found:
            #    posted.append( (out_type,rep_tech,fid) )
            
            fileDict = dxencode.description_from_fid(fid,properties=True)
            #file_name = dxencode.file_path_from_fid(fid,projectToo=False)
            acc_key = dxencode.dx_property_accesion_key('https://www.encodeproject.org') # prefer production accession
            # check file properties
            accession = ''
            if "properties" in fileDict:
                if acc_key not in fileDict["properties"] and self.server_key != 'www':   # test, beta replicated from www
                    acc_key = dxencode.dx_property_accesion_key(self.server)
                if acc_key in fileDict["properties"]:
                    accession = fileDict["properties"][acc_key]
                    if verbose:
                        print >> sys.stderr, "* DEBUG   Accession: " + accession
                    f_obj =  self.enc_file_find_by_dxid(fid)  # Look by alias first incase ENCFF v. TSTFF mismatch
                    if f_obj != None: # Verifiably posted
                        if f_obj.get('status') == 'revoked':
                            if verbose:
                                print >> sys.stderr, "* DEBUG   Found REVOKED by alias"
                            self.revoked.append(accession)
                            self.enc_file_add_alias(fid,accession,f_obj,remove=True,test=test)
                            ret = dxencode.dx_file_set_property(fid,'revoked_accession',accession,test=test)
                            if not test and ret == accession:
                                dxencode.dx_file_set_property(fid,'accession',null,test=test)
                                print "  * Marked DX property with revoked_accession."
                        else:
                            self.found[fid] = f_obj
                            if f_obj.get('status') != 'upload failed': 
                                posted.append( (out_type,rep_tech,fid) )
                            # TODO: Compare accessions.  But not ENCFF to TSTFF
                            if verbose:
                                print >> sys.stderr, "* DEBUG   Found by alias"
                                if f_obj.get('status') == 'upload failed': 
                                    print >> sys.stderr, "* DEBUG   But needs to be reposted"
                        continue

                    # look by accession
                    if verbose:
                        print >> sys.stderr, "* DEBUG   Not found by alias: dnanexus:" + fid
                    f_obj = self.enc_lookup_json( 'files/' + accession,must_find=False)
                    if f_obj != None: # Verifyably posted
                        if f_obj.get('status') == 'revoked':
                            if verbose:
                                print >> sys.stderr, "* DEBUG   Found REVOKED by accession"
                            self.revoked.append(accession)
                            ret = dxencode.dx_file_set_property(fid,'revoked_accession',accession,test=test)
                            if not test and ret == accession:
                                dxencode.dx_file_set_property(fid,'accession',null,test=test)
                                print "  * Marked DX property with revoked_accession."
                        else:
                            self.found[fid] = f_obj
                            if f_obj.get('status') != 'upload failed': 
                                posted.append( (out_type,rep_tech,fid) )
                            self.enc_file_add_alias(fid,accession,f_obj,test=test)
                            if verbose:
                                print >> sys.stderr, "* DEBUG   Found by accession: " + accession
                                if f_obj.get('status') == 'upload failed': 
                                    print >> sys.stderr, "* DEBUG   But needs to be reposted"
                        continue
                    else:
                        if verbose:
                            print >> sys.stderr, "* DEBUG   Not found by accession: " + accession
                        
            if accession == '':
                if verbose:
                    print "* DEBUG Accession not found."
                # No accession in properties, but try to match by fid anyway.
                #f_obj =  self.find_in_encode(fid,verbose)
                f_obj =  self.enc_file_find_by_dxid(fid)
                if f_obj != None:
                    accession = f_obj['accession']
                    if f_obj.get('status') == 'revoked':
                        if verbose:
                            print >> sys.stderr, "* DEBUG   Found REVOKED by alias"
                        self.revoked.append(accession)
                        self.enc_file_add_alias(fid,accession,f_obj,remove=True,test=test)
                    else:
                        self.found[fid] = f_obj
                        if f_obj.get('status') != 'upload failed': 
                            posted.append( (out_type,rep_tech,fid) )
                        # update dx file property
                        if accession.startswith('ENCFF'):  # Only bother if it is the real thing.  Ambiguity with 'TSTFF'
                            ret = dxencode.dx_file_set_property(fid,'accession',accession,add_only=True,test=test)
                            if not test:
                                if ret == accession:
                                    print "  * Updated DX property with accession."
                        if verbose:
                            print >> sys.stderr, "* DEBUG   Found by alias"
                            if f_obj.get('status') == 'upload failed': 
                                print >> sys.stderr, "* DEBUG   But needs to be reposted"
                    continue

                if verbose:
                    print >> sys.stderr, "* DEBUG   Not found by alias: dnanexus:" + fid

            # Check the way back machine of file name and size.
            if self.way_back_machine:
                f_obj =  self.find_files_using_way_back_machine(fid,self.exp_files) 
                if f_obj != None:
                    self.found[fid] = f_obj
                    if f_obj.get('status') != 'upload failed': 
                        posted.append( (out_type,rep_tech,fid) )
                    # Update ENC and DX:
                    accession = f_obj['accession']
                    if accession.startswith('ENCFF'):  # Only bother if it is the real thing.  Ambiguity with 'TSTFF'
                        ret = dxencode.dx_file_set_property(fid,'accession',accession,add_only=True,test=test)
                        if not test:
                            if ret == accession:
                                print "  * Updated DX property with accession."
                        self.enc_file_add_alias(fid,accession,f_obj,test=test)
                                
                    if verbose:
                        print >> sys.stderr, "* DEBUG   Found using 'way back machine'."
                        if f_obj.get('status') == 'upload failed': 
                            print >> sys.stderr, "* DEBUG   But needs to be reposted"
                    continue
                if verbose:
                    print >> sys.stderr, "* DEBUG   Not Found using 'way back machine'.  VERBOSE..."
                    f_obj =  self.find_files_using_way_back_machine(fid,self.exp_files,verbose=True)

        if verbose:
            print >> sys.stderr, "Posted files:"
            print >> sys.stderr, json.dumps(posted,indent=4)
        return posted
        
    def find_needed_files(self,files_expected,test=True,verbose=False):
        '''Returns the tuple list of files that NEED to be posted to ENCODE.'''
        needed = []
        posted = self.find_posted_files(files_expected,test=test,verbose=verbose)
        #verbose = True
        # FIXME: special case to get around already updated aliases
        #self.revoked.extend(['ENCFF905HBB','ENCFF374TIF','ENCFF659SBG','ENCFF308DBJ','ENCFF840CVB'])
        # FIXME: special case to get around already updated aliases
        
        # Note that find_posted_files() fills in self.found
        for (out_type, rep_tech, fid, QC_only) in files_expected:
            if not QC_only and (out_type, rep_tech, fid) not in posted:
                needed.append( (out_type,rep_tech,fid, False) )
                if verbose:
                    print >> sys.stderr, "* DEBUG post needed for:" + fid
                continue
            # Figure out if a patch is needed
            # WARNING: This will only work if the revoked files have the old dx aliases, and this run will update those...
            #          so patching derived from when earlier files were revoke will only work on the first non-test run.
            if not QC_only and (out_type, rep_tech, fid) in posted:
                patch_needed = False
                if fid in self.found.keys():
                    f_obj = self.found[fid]
                    derived_from = f_obj.get('derived_from')
                    if derived_from != None and isinstance(derived_from,list):
                        for file_url in derived_from:
                            revoked_acc = file_url.split('/')[2]
                            if revoked_acc in self.revoked:
                                patch_needed = True
                                break
                if patch_needed:
                    needed.append( (out_type,rep_tech,fid, False) )
                    if verbose:
                        print >> sys.stderr, "* DEBUG patch of derived_from is needed for:" + fid
                    continue
            # Figure out if QC was already posted
            qc_json = self.dx_qc_json_get(fid)
            if not qc_json:
                continue
            for qc_key in qc_json.keys():
                if qc_key not in self.QC_SUPPORTED:
                    continue
                qc_metric = self.enc_qc_metric_find(fid,qc_key)
                if not qc_metric or qc_metric == None:
                    if verbose:
                        print >> sys.stderr, "* DEBUG QC_only is needed for:" + fid
                    needed.append( (out_type,rep_tech,fid, True) )
                    break
                else:
                    f_obj = self.found.get(fid)
                    if f_obj != None and "accession" in f_obj:
                        file_ref = "/files/%s/" %f_obj["accession"]
                        if file_ref not in qc_metric["quality_metric_of"]:
                            if verbose:
                                print >> sys.stderr, "* DEBUG QC_only patch is needed for:" + fid
                            needed.append( (out_type,rep_tech,fid, True) )
                            break
        if verbose:
            print >> sys.stderr, "Needed files:"
            print >> sys.stderr, json.dumps(needed,indent=4)
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
            print >> sys.stderr, "sw_versions:"
            print >> sys.stderr, json.dumps(sw_versions,indent=4)
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
            print >> sys.stderr, "app_version:"
            print >> sys.stderr, json.dumps(app_version,indent=4)
        return app_version
        
    def input_exception(self,inp_fid):
        '''Returns True if this is one of a limit number of input files we do not track in encodeD.'''
        # TODO: move specifics to json at top of file.
        if self.exp_type == "dna-me" and dxencode.file_path_from_fid(inp_fid).endswith("_techrep_bismark_map_report.txt"):
            #print "  Ignoring file: " + dxencode.file_path_from_fid(inp_fid)
            return True
        return False
                
    def find_derived_from(self,fid,job,verbose=False):
        '''Returns list of accessions a file is drived from based upon job inouts.'''
        input_accessions = []
        input_file_count = 0
        #verbose=True  # NOTE: verbose can help find the missing accessions  # excessive verbosity helped debugging recovery.py
        # NOTE: What to do about concatenated files (that went through 'concat_fastqs')?
        #       Can punt for now since recovery can append to derived_from
        
        # Find all expected inputs from the job
        file_inputs = []
        for inp in job["input"].values():
            if type(inp) not in [ dict, list ]:
                continue # not a file input
            if type(inp) == dict:
                file_inputs.append(inp)
                continue
            for item in inp:
                if type(item) == dict:
                    file_inputs.append(item)
        if verbose:
            print >> sys.stderr, "* derived from: Expecting %d input files." % len(file_inputs)
            
        # For each file input, verify it is for a file and then look for an accession.
        for inp in file_inputs:
            if not type(inp) == dict:
                print type(inp)
                continue # not a file input
            dxlink = inp.get("$dnanexus_link")
            if dxlink == None:
                continue
            if not type(dxlink) == dict:
                inp_fid = dxlink
            else:
                inp_fid = dxlink.get("id")
            if self.input_exception(inp_fid):  # Look for exceptions
                continue
            input_file_count += 1
            try:
                inp_obj = dxencode.description_from_fid(inp_fid,properties=True)
            except:
                print "WARNING: can't find "+ fid # may try to append derived_from below.
                continue
            if verbose: 
                print >> sys.stderr, "* derived from: " + inp_fid + " " + inp_obj["project"] + ":" + \
                                                                        dxencode.file_path_from_fid(inp_fid)
            if inp_obj == None:
                continue
                
            # First best place to look for accession: properties
            self.genome = self.find_genome_annotation(fid,inp_obj)
            acc_key = dxencode.dx_property_accesion_key(dxencode.PRODUCTION_SERVER) # Always prefer production accession
            accession = None
            if "properties" in inp_obj:
                if verbose:
                    print >> sys.stderr, "Input file properties... looking for '"+acc_key+"'" 
                    print >> sys.stderr, json.dumps(inp_obj["properties"],indent=4)
                # older runs don't have accession in properties
                if acc_key in inp_obj["properties"]:
                    accession = inp_obj["properties"][acc_key]
                    if verbose:
                        print >> sys.stderr, "Found accession."
                    # Must check if file exists!!
                    file_obj = self.enc_lookup_json( 'files/' + accession,must_find=False)
                    if file_obj == None:
                        if verbose:
                            print >> sys.stderr, "Accession found but file not on '"+self.server_key+"'"
                        accession = None
                    else:
                        input_accessions.append(accession)                        
                            
                # may need to look for server specific accession
                if accession == None and self.server_key != dxencode.PRODUCTION_SERVER:
                    acc_key = dxencode.dx_property_accesion_key(self.server)  # demote to server's accession
                    if verbose:
                        print >> sys.stderr, "Now looking for '"+acc_key+"' on '"+self.server_key+"'" 
                    if acc_key in inp_obj["properties"]:
                        accession = inp_obj["properties"][acc_key]
                        if verbose:
                            print >> sys.stderr, "Found accession for '"+self.server_key+"'."
                        # Must check if file exists!!
                        file_obj = self.enc_lookup_json( 'files/' + accession,must_find=False)
                        if file_obj == None:
                            if verbose:
                                print >> sys.stderr, "Accession found but file not on '"+self.server_key+"'"
                            accession = None
                        else:
                            input_accessions.append(accession)                        
            
            # Now start the tricks to in the vain hope to turn up an accession
            if accession == None:
                if verbose:
                    print >> sys.stderr, "Accession not found in properties or not in ENCODEd." 
                parts = inp_obj["name"].split('.')
                ext = parts[-1]
                if ext in ["gz","tgz"]:
                    ext = parts[-2]
                if ext in self.REFERENCE_EXTENSIONS and inp_obj["project"] != self.REFERERNCE_PROJECT_ID:  
                    # FIXME: if regular file extension is the same as ref extension, then trouble!
                    print "Reference file (ext: '"+ext+"') in non-ref project: " + inp_obj["project"]
                    sys.exit(0)

                # if file name is primary input (fastq) and is named as an accession
                if inp_obj["name"].startswith("ENCFF"): # Not test version 'TSTFF'!
                    if ext in self.PRIMARY_INPUT_EXTENSION or self.test:
                        root = parts[0]
                        for acc in root.split('_'): #usually only one
                            if acc.startswith("ENCFF"):
                                if len(acc) == 11:
                                    accession = acc
                                    input_accessions.append(accession)
                                elif '-' in acc:
                                    for acc_part in acc.split('-'):
                                        if acc_part.startswith("ENCFF") and len(acc_part) == 11:
                                            accession = acc_part
                                            input_accessions.append(accession)
                                            
                # no accession but test run of splashdown and not fastq just assume a real run would have posted by now.
                if accession == None and self.test and self.TOOL_IS == 'splashdown' and \
                    ext not in self.PRIMARY_INPUT_EXTENSION:
                    # since we are testing splashdown, lets assume a real run would have posted the previous files
                    if verbose:
                        print >> sys.stderr, "Test run so using 'ENCFF00FAKE' as accession of derived from file."
                    accession = "ENCFF00FAKE"
                    input_accessions.append(accession)
                        
            # Still no accession, try the way back machine
            if accession == None and self.way_back_machine and self.exp_files != None:
                if verbose:
                    print >> sys.stderr, "Accession still not found.  Trying the 'way back machine'." 
                f_obj =  self.find_files_using_way_back_machine(inp_fid,self.exp_files)
                if f_obj != None and "accession" in f_obj:
                    accession = f_obj["accession"]
                    input_accessions.append(accession)
    
        # All file_inputs were tested and we hope that the number of accessions match the number of input files.
        if len(input_accessions) < input_file_count:
            if self.test:
                if self.way_back_machine and self.TOOL_IS == 'recovery':
                    print "WARNING: not all input files are accounted for in 'derived_from' so just appending! Found %d of %d." % \
                                                                        (len(input_accessions),input_file_count)
                    input_accessions.append(self.APPEND_FLAG)
                else:
                    print "WARNING: not all input files are accounted for in 'derived_from'! Found %d of %d." % \
                                                                        (len(input_accessions),input_file_count)
                    print json.dumps(input_accessions,indent=4)
            else:
                if self.way_back_machine and self.TOOL_IS == 'recovery':
                    print "WARNING: not all input files are accounted for in 'derived_from' so just appending!"
                    input_accessions.append(self.APPEND_FLAG)
                else:
                    print "ERROR: not all input files are accounted for in 'derived_from'! Found %d of %d." % \
                                                                        (len(input_accessions),input_file_count)
                    print json.dumps(input_accessions,indent=4)
                    sys.exit(1)
        
        # Now that we have the full derived_from, we can remove some         
        # UGLY special cases:
        derived_from = []
        for acc in input_accessions:
            if not acc.endswith("chrom.sizes"):  # Chrom.sizes is just a bit too much info
                derived_from.append(acc)
        
        if verbose:
            print >> sys.stderr, "Derived files: for " + dxencode.file_path_from_fid(fid)
            print >> sys.stderr, json.dumps(derived_from,indent=4)
        return derived_from


    def add_encoded_info(self,obj,rep_tech,fid,verbose=False):
        '''Updates an object with information from encoded database.'''
        obj['lab'] = '/labs/encode-processing-pipeline/' # self.exp['lab']['@id'] Now hard-coded
        obj['award'] = '/awards/U41HG006992/'  # self.exp['award']['@id']
        #obj['submitted_by'] = '/users/8b1f8780-b5d6-4fb7-a5a2-ddcec9054288/' # or 'ENCODE Consortium' or '/labs/encode-consortium/' 

        # Find replicate info
        #if rep_tech.startswith("rep") and len(rep_tech) == 6:
        (br,tr) = (None,None)
        if rep_tech.startswith("reps"):
            br_tr = rep_tech[4:]
            (br,tr) = br_tr.split('_')
            tr = tr.split('.')[-1]       # TODO: Get buy in for bio_rep files being associated with the last tech_rep.
        elif rep_tech.startswith("rep"):
            br_tr = rep_tech[3:]
            (br,tr) = br_tr.split('_')
        if br != None and tr != None:
            full_mapping = dxencode.get_full_mapping(self.exp_id,self.exp)
            mapping = dxencode.get_replicate_mapping(self.exp_id,int(br),int(tr),full_mapping)
            obj['replicate'] = mapping['replicate_id']

        if verbose:
            print >> sys.stderr, "After adding encoded info:"
            print >> sys.stderr, json.dumps(obj,indent=4)
        return obj


    def dx_qc_json_get(self,fid,verbose=False):
        '''Get the qc metrics blob from dx details, with very special exceptions'''
        qc_json = {}
        descr = dxencode.description_from_fid(fid)
        file_name = descr['name']
        # NOTE: So far only special case for getting qc_from_file 
        if not self.qc_from_file or (not file_name.endswith('_star_genome.bam') and not file_name.endswith('_star_anno.bam')):
            qc_json = dxencode.dx_file_get_details(fid)
            if qc_json and "QC" in qc_json:  # Not likely but QC json could be subsection of details
                qc_json = qc_json["QC"]
            
        # Try to make a qc_json in a very special case
        if not qc_json:
            if not file_name.endswith('_star_genome.bam') and not file_name.endswith('_star_anno.bam'):
                return {}
            folder = descr['folder']
            qc_fid = dxencode.find_file(folder + "/*_star_Log.final.out",self.proj_id,multiple=False,recurse=False)
            if qc_fid == None:
                return {}
                
            # 1) copy locally to tmp.txt  (should be fast)
            if not os.path.isfile("tmp/"+qc_fid+".txt") :
                try:
                    os.system('mkdir -p tmp > /dev/null 2>&1')
                    dxpy.download_dxfile(qc_fid, "tmp/"+qc_fid+".txt")  # should be short
                except:
                    print >> sys.stderr, "ERROR: Unable to download '"+folder + "/*_star_Log.final.out'."
                    return {}

            # 2) use qc_metrics.py to shlurp in the json 
            qc_parser = "~/tim/long-rna-seq-pipeline/dnanexus/tools/qc_metrics.py"
            err, out = commands.getstatusoutput(qc_parser + " -n STAR_log_final --json -f tmp/"+qc_fid+".txt 2> /dev/null")
            # ignore stderr as it is used to echo results
            if len(out) > 0:
                qc_json = { "STAR_log_final": json.loads(out) } 
            
        if verbose and qc_json:
            print >> sys.stderr, "DX 'qc_json':"
            print >> sys.stderr, json.dumps(qc_json,indent=4,sort_keys=True)
        return qc_json

    def enc_qc_metric_find(self,fid,qc_key,job_id=None,collection=None,must_find=False):
        '''Returns the qc_mtrics object from either cache or encodeD.'''
        qc_alias = self.qc_metric_make_alias(fid,qc_key,job_id) 
        if "exp" in self.obj_cache and qc_alias in self.obj_cache["exp"]:
            return self.obj_cache["exp"][qc_alias]

        # What collection (schema type) will this metric belong to?
        if collection == None:
            collection = self.qc_metric_schema_type(qc_key)

        qc_metric = self.enc_lookup_json(qc_alias,must_find=must_find)
        if qc_metric != None:
            self.obj_cache["exp"][qc_alias] = qc_metric
        return qc_metric
        

    def qc_props_fill(self,dx_qc_obj,qc_faq):
        '''Fills in the qc properties from a dx_qc_obj, taking into account any inclues, excludes, etc.'''
        enc_qc_props = {}
        if "singleton" in qc_faq:
            enc_qc_props[qc_faq["singleton"]] = dx_qc_obj
        elif "props" in qc_faq:
            for key in qc_faq['props'].keys():
                enc_qc_props[qc_faq['props'][key]] = dx_qc_obj[key]
        elif "include" in qc_faq:
            for key in qc_faq['include']:
                enc_qc_props[key] = dx_qc_obj[key]
        elif "exclude" in qc_faq:
            for key in dx_qc_obj.keys():
                if key not in qc_faq['exclude']:
                    enc_qc_props[key] = dx_qc_obj[key]
        else:
            enc_qc_props = dx_qc_obj # All props in DX object shall be placed in enc object
        return enc_qc_props


    def qc_metric_attachment(self,fid,blob_def):
        '''Returns a QC metrics attachment property for posting an attachment with the QC object.'''
        attachment = {}
        blob_fid = None
        # "blob": { "pattern": "/*_mad_plot.png" }, 
        # from base64 import b64encode

        # Need to find the dx file (fid)
        # should be in same folder as fid
        descr = dxencode.description_from_fid(fid)
        
        path = dxencode.file_path_from_fid(fid)
        folder = descr['folder']
        file_name = descr['name']
        if file_name.endswith(blob_def['pattern'][2:]):  # is this the file?
            blob_fid = fid
            blob_name = file_name
            blob_size = descr['size']
        if blob_fid == None:
            blob_fid = dxencode.find_file(folder + blob_def['pattern'],self.proj_id,multiple=False,recurse=False)
            if blob_fid == None:
                print >> sys.stderr, "ERROR: Failed to find blob attachment file: " +folder + blob_def['pattern']
                sys.exit(1) # Make this a halt to see what is going on
            blob_descr = dxencode.description_from_fid(blob_fid)
            blob_name = blob_descr['name']
            blob_size = blob_descr['size']
        
        # What type of blob is this?
        if "mime_ext" in blob_def:
            blob_name = blob_name + "." + blob_def["mime_ext"]
        blob_mime = self.file_mime_type(blob_name)
        if blob_mime == None:
            blob_mime = "text/plain"  # Can we just do this?
            print >> sys.stderr, "ERROR: Unrecognizible mime_type for file: " + blob_name
            sys.exit(1) 

        # Stream the file and fill in attachment
        try:
            with dxpy.open_dxfile(blob_fid) as stream:
                attachment = {
                    'download': blob_name, #Just echoes the given filename as the download name
                    'size':     blob_size,
                    'type':     blob_mime,
                    'href': 'data:%s;base64,%s' % (blob_mime, b64encode(stream.read()))
                }
            #print "  * Created attachemnt '"+blob_name+"'."
            return attachment
        except:
            print >> sys.stderr, "ERROR: Unable to open and read '"+blob_name+"' as stream."
        return None


    def qc_metric_make_alias(self,fid,qc_key,job_id=None):
        '''Returns a qc_metric alias'''
        if job_id == None:
            try:
                file_dict = dxencode.description_from_fid(fid)
                job_id = file_dict["createdBy"]["job"]
            except:
                print >> sys.stderr, "ERROR: could not find job for "+fid
                return None
        qc_type = qc_key
        # The alias must have qc_key, NOT "type" because one job could have 2 different samtools_flagstats!
        #if qc_key in self.QC_SUPPORTED and "type" in self.QC_SUPPORTED[qc_key]:
        #    qc_type = self.QC_SUPPORTED[qc_key]["type"]
        return "dnanexus:qc."+qc_type+'.'+job_id
        
        
    def qc_metric_schema_type(self,qc_key):
        '''Returns the schema type (collection) a qc_metric belongs to'''
        qc_suffix = "_quality_metric"
        if qc_key in self.QC_SUPPORTED and "type" in self.QC_SUPPORTED[qc_key]:
            return self.QC_SUPPORTED[qc_key]["type"] + qc_suffix
        else: 
            return qc_key.lower() + qc_suffix


    def qc_metric_files(self,qc_faq,fid,verbose=False):
        '''Returns a list of file accessions that this qc metric is supposed to relate to.'''
        if "files" not in qc_faq:
            return []
        qc_fids = []
        if "results" in qc_faq["files"].keys():  # the fid could be an attachment or a posted file
            qc_fids.append(fid) # if fid is for the attachment, the conversion to an accession will eliminate it.
        if "inputs" in qc_faq["files"].keys():
            inp_names = qc_faq["files"]["inputs"]
            job = dxencode.job_from_fid(fid)
            job_inputs = job.get('input')
            for name in inp_names:
                if name in job_inputs:
                    dx_inp = job_inputs[name]
                    if verbose:
                        print >> sys.stderr, "Job inputs:"
                        print >> sys.stderr, json.dumps(dx_inp,indent=4,sort_keys=True)
                    if isinstance(dx_inp,list):
                        for link in dx_inp:
                            qc_fids.append(link["$dnanexus_link"])
                    else:
                        qc_fids.append(dx_inp["$dnanexus_link"])
        qc_accs = []
        for qc_fid in qc_fids:
            acc = self.file_get_accession(qc_fid,fake_if_needed=self.test)
            if acc != None:
                qc_accs.append('/files/'+acc+'/')
        if verbose:
            print >> sys.stderr, "Expect qc_metric['quality_metric_of'] fids:"
            print >> sys.stderr, json.dumps(qc_fids,indent=4,sort_keys=True)
            print >> sys.stderr, "Expect qc_metric['quality_metric_of']:"
            print >> sys.stderr, json.dumps(qc_accs,indent=4,sort_keys=True)

        return qc_accs


    def enc_qc_metric_find_or_create(self,qc_key,qc_obj,qc_faq,fid,job_id,step_run_id,test=True,verbose=False):
        '''Finds or creates the 'qc_metric' encoded object that will be attached to the step_run and/or file.'''
        #verbose=True
        qc_metric = None
        qc_patch = {}
        
        qc_props = self.qc_props_fill(qc_obj,qc_faq)
        if not qc_props:
            return None

        collection = self.qc_metric_schema_type(qc_key)
        qc_files = self.qc_metric_files(qc_faq,fid)
        qc_alias = self.qc_metric_make_alias(fid,qc_key,job_id)        
        
        qc_metric = self.enc_qc_metric_find(fid,qc_key,job_id,collection)
        if qc_metric != None:
            # How to tell if it was just created???
            # The step_run will only be an alias on just created objects
            if qc_metric['step_run'] == step_run_id:
                print "  - Found qc_metric: '%s'" % qc_alias
                self.obj_cache["exp"][qc_alias] = qc_metric
                for prop in qc_props.keys():
                    if None == qc_metric.get(prop):
                        print >> sys.stderr, "ERROR: Expecting '"+prop+"' in qc_metric."
                        qc_patch=qc_props
                        #break
                        sys.exit(1) # Until a legit case comes along, better exit and figure it out
                    elif qc_metric[prop] != qc_props[prop]:
                        print >> sys.stderr, "ERROR: qc_metric['"+prop+"'] expecing <"+str(qc_props[prop])+ \
                                                                                    ">, but found <"+str(qc_metric[prop])+">."
                        qc_patch=qc_props
                        #break
                        sys.exit(1) # Until a legit case comes along, better exit and figure it out
                
        if qc_metric:
            # update files in quality_metric_of list
            # Note ignoring patch of attachment at this time
            if len(qc_files) > 0:
                for qc_file in qc_files:
                    if qc_file not in qc_metric['quality_metric_of']:
                        qc_metric['quality_metric_of'].append(qc_file)
                        qc_patch['quality_metric_of'] = qc_metric['quality_metric_of']                        
            
            # Now patch if necessary        
            if qc_patch:
                if verbose:
                    print >> sys.stderr, "ENCODEd 'qc_metric' patch:"
                    print >> sys.stderr, json.dumps(qc_patch,indent=4,sort_keys=True)
                    
                self.obj_cache["exp"][qc_alias] = qc_metric
                if test:
                    print "  * Would patch qc_metric: '%s'" % qc_alias
                else:
                    try:
                        patched_obj = dxencode.encoded_patch_obj(qc_alias,qc_patch, self.server, self.authid, self.authpw)
                    except:
                        print "Failed to patch qc_metric: '%s'" % qc_alias
                        sys.exit(1)
                    print "  * Patched qc_metric: '%s'" % qc_alias
                    #self.qc_metric_created += 1
                
        # No found metric so create
        if qc_metric == None:
            qc_metric = qc_props
                
            if len(qc_files) > 0:
                qc_metric['quality_metric_of'] = qc_files
            qc_metric['step_run'] =  step_run_id
            qc_metric['assay_term_name'] = self.exp['assay_term_name']
            qc_metric['assay_term_id']   = self.exp['assay_term_id'] # self.exp_type 
            qc_metric['aliases'] = [ qc_alias ]
            if verbose:
                print >> sys.stderr, "ENCODEd 'qc_metric' before attachment:"
                print >> sys.stderr, json.dumps(qc_metric,indent=4,sort_keys=True)

            blob_msg = ""
            if "blob" in qc_faq:
                attachment = self.qc_metric_attachment(fid,qc_faq["blob"])
                if attachment != None:
                    qc_metric['attachment'] = attachment
                    blob_msg = " (with attachment)"
                
            self.obj_cache["exp"][qc_alias] = qc_metric
            if test:
                print "  * Would post qc_metric%s: '%s'" % (blob_msg,qc_alias)
                #print json.dumps(qc_metric,indent=4,sort_keys=True)
            else:
                try:
                    posted_obj = dxencode.encoded_post_obj(collection,qc_metric, self.server, self.authid, self.authpw)
                except:
                    print "Failed to post qc_metric%s: '%s'" % (blob_msg,qc_alias)
                    sys.exit(1)
                print "  * Posted qc_metric%s: '%s'" % (blob_msg,qc_alias)
                #self.qc_metric_created += 1
                
        return qc_metric
        
    def find_qc_key(self,key,payload):
        '''QC instructions may differ based upon file type.'''
        for qc_key in self.QC_SUPPORTED.keys():
            qc_faq = self.QC_SUPPORTED[qc_key]
            if "only_for" in qc_faq:
                if qc_faq.get("type") != key:
                    continue
                if isinstance(qc_faq["only_for"],str):
                    if payload["submitted_file_name"].endswith(qc_faq["only_for"]):
                        return qc_key
                elif isinstance(qc_faq["only_for"],list):
                    for one_ending in qc_faq["only_for"]:
                        if payload["submitted_file_name"].endswith(one_ending):
                            return qc_key
        return key
            
        
    def handle_qc_metrics(self,fid,payload,test=True,verbose=False):
        '''After posting a file, post or patch any qc_metric objects associated with it.'''
        # Strategy:
        # 1) ALL qc metrics MUST be associated with a file - file may not be posted ('QC_only')
        # 2) read file's "details" json string
        # 3) For each key in qc_json object: if it matches known type: pass to specialized routine to
        #    a) Find existing by dnanexus:qc.qc_key.job-2342345:
        #    b) If found, patch as necessary (expect to need to add file to files list)
        #    c) If not found, create, optionally attach blob, and post it
        #verbose=True
        
        # Needed identifiers
        step_run_id = payload["step_run"]
        step_run_alias = step_run_id.split('/')[-1]
        if step_run_alias.startswith('dnanexus:'):
            job_id = step_run_alias.split(':')[-1]
        else:
            job = dxencode.job_from_fid(fid)
            job_id = job.get('id')
            #step_run_alias = 'dnanexus:' + job_id

        # if file level qc_metric exist, then use those
        qc_for_file = {}
        qc_json = self.dx_qc_json_get(fid)
        if not qc_json:
            return None
        for key in qc_json:
            qc_key = self.find_qc_key(key,payload)
            if verbose:
                print >> sys.stderr, "  * Working on qc_key: '%s'" % (qc_key)                
            if qc_key not in self.QC_SUPPORTED:
                continue
            qc_faq = self.QC_SUPPORTED[qc_key]
            qc_metric = self.enc_qc_metric_find_or_create(qc_key,qc_json[key],qc_faq,fid,job_id,step_run_id, \
                                                                                            test=test,verbose=verbose)
            if qc_metric != None:
                qc_for_file[key] = qc_json[key]
                
        if qc_for_file and verbose:
            print >> sys.stderr, "qc_for_file:"
            print >> sys.stderr, json.dumps(qc_for_file,indent=4,sort_keys=True)
            
        return len(qc_for_file)

    def enc_step_version_find(self,step_ver_alias,verbose=False):
        '''Finds the 'analysis_step_version' encoded object used in creating the file.'''
        step_ver = None
        
        if step_ver_alias in self.obj_cache:
            return self.obj_cache[step_ver_alias]
        step_ver = self.enc_lookup_json( 'analysis-step-versions/' + step_ver_alias,must_find=True)
        if step_ver:
            self.obj_cache[step_ver_alias] = step_ver
            print "  - Found step_ver: '%s'" % step_ver_alias
            if verbose:
                print json.dumps(step_ver,indent=4,sort_keys=True)
        return step_ver


    def enc_step_run_find_or_create(self,job,dxFile,rep_tech,test=False,verbose=False):
        '''Finds or creates the 'analysis_step_run' encoded object that actually created the file.'''
        #verbose=True
        
        step_run = None
        job_id = job.get('id')
        step_alias = 'dnanexus:' + job_id
        if "exp" in self.obj_cache and step_alias in self.obj_cache["exp"]:
            step_run = self.obj_cache["exp"][step_alias]
        else:
            step_run = self.enc_lookup_json( 'analysis-step-runs/' + step_alias,must_find=False)
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
            
            # get applet and version
            dx_app_id = job.get('applet')
            dx_app_ver = job.get('step_parent_app_version') # Could have be different from child because this is a step parent.
            if dx_app_ver == None:
                dx_app_ver = self.find_app_version(dxFile)
            if dx_app_ver and 'version' in dx_app_ver:
                dx_app_ver = str( dx_app_ver.get('version') )
                if dx_app_ver[0] == 'v':
                    dx_app_ver = dx_app_ver[1:]
            # FIXME: UGLY special case
            if dx_app_ver.startswith('0.'):  # If posting results generated while still developing, assume v1.0.0
                dx_app_ver = "1.0.0"
            # FIXME: UGLY temporary special case!!!
            if dx_app_name  == "rampage-peaks" and dx_app_ver == "1.0.1":
                dx_app_ver = "1.1.0"  # Because the version was supposed to bump the second digit if a tool changes.
            # FIXME: UGLY temporary special case!!!
            if not dx_app_ver or not isinstance(dx_app_ver, str) or len(dx_app_ver) == 0:
                print "ERROR: cannot find applet version %s in the log" % ( type(dx_app_ver) )
                sys.exit(0)
                
            # get analysis_step_version (aka step_ver):
            # NOTE: the special dnanexus: alias.  Because dx_ids will differ between projects, and the identical app can be
            #       rebuilt, the dnanexus: alias will instead be the app_name and first 2 digits of the app_version!
            step_ver_alias = "dnanexus:" + dx_app_name + '-v-' + '-'.join(dx_app_ver.split('.')[:-1])
            step_ver = self.enc_step_version_find(step_ver_alias)
            step_run['analysis_step_version'] = '/analysis-step-versions/' +step_ver_alias
            
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
            if verbose:
                print "DX JOB: originalInputs:"
                print json.dumps(inputs,indent=4,sort_keys=True)
            for name in inputs.keys():
                #if isinstance(inputs[name],str) or isinstance(inputs[name],int) or isinstance(inputs[name],unicode):
                if not isinstance(inputs[name],dict):  # not class file? and not class:array:file 
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
            # analysis_step? No, can get that from analysis_step_version
            # pipeline?      No, can get that from analysis_step
            notes["notes_version"] = "2"
            notes["dx_app_name"] = dx_app_name
            notes['dx_app_version'] = dx_app_ver
            notes["dx_app_id"] = dx_app_id
            notes["step_name"] = step_ver['analysis_step']
            notes["dx_analysis_id"] = job.get('analysis')
            if "ana_id" not in self.obj_cache["exp"]:
                self.obj_cache["exp"]["ana_id"] = [ notes["dx_analysis_id"] ]
            elif notes["dx_analysis_id"] not in self.obj_cache["exp"]["ana_id"]:
                self.obj_cache["exp"]["ana_id"].append( notes["dx_analysis_id"] )
            notes["dx_project_id"] = self.proj_id
            notes["dx_project_name"] = self.proj_name
            notes["dx_cost"] = "$" + str(round(job['totalPrice'],2))
            duration = dxencode.format_duration(job.get('startedRunning')/1000,job.get('stoppedRunning')/1000)
            notes["duration"] = duration
            step_run["notes"] = json.dumps(notes)
            
            # Now post this new object
            self.obj_cache["exp"][step_alias] = step_run
            if test:
                print "  * Would post step_run: '%s'" % step_alias
            else:
                try:
                    #step_run["@type"] = ["item", "analysis_step_run"]
                    posted_step_run = dxencode.encoded_post_obj('analysis_step_run',step_run, self.server, self.authid, self.authpw)
                except:
                    print "Failed to post step_run: '%s'" % step_alias
                    sys.exit(1)
                print "  * Posted step_run: '%s'" % step_alias
                self.step_runs_created += 1

        if step_run and verbose:
            print >> sys.stderr, "ENCODEd 'analysis_step_run':"
            print >> sys.stderr, json.dumps(step_run,indent=4)
            if "notes" in step_run:
                step_run_notes = json.loads(step_run.get("notes"))
                print >> sys.stderr, "ENCODEd 'analysis_step_run[notes]':"
                print >> sys.stderr, json.dumps(step_run_notes,indent=4)
        return step_run
        
        
    def find_step_parent(self,fid,child_job,step_child,step_parent,verbose=False):
        '''Returns the parent job for a child that is not recoreded to encodeD.'''
        #print >> sys.stderr, "Step-child %s points to %s" % (step_child,step_parent)
        #print >> sys.stderr, json.dumps(child_job,indent=4,sort_keys=True)
        #sys.exit(1)
        
        # Look for parent by walking up derived from files
        file_inputs = []
        for inp in child_job["input"].values():
            if type(inp) not in [ dict, list ]:
                continue # not a file input
            if type(inp) == dict:
                file_inputs.append(inp)
                continue
            for item in inp:
                if type(item) == dict:
                    file_inputs.append(item)
            
        # For each file input, verify it is for a file and then look for an accession.
        parent_job = None
        parent_fid = None
        for inp in file_inputs:
            if not type(inp) == dict:
                print type(inp)
                continue # not a file input
            dxlink = inp.get("$dnanexus_link")
            if dxlink == None:
                continue
            if not type(dxlink) == dict:
                inp_fid = dxlink
            else:
                inp_fid = dxlink.get("id")
            try:
                job = dxencode.job_from_fid(inp_fid)
            except:
                print "WARNING: can't find parent_job for "+ fid # may try to append derived_from below.
                continue
            if job.get('executableName') == step_parent:
                parent_job = job
                parent_fid = inp_fid
                break
        if parent_job == None:
            print "ERROR: Step-child %s cannot find it's step-parent %s" % (step_child,step_parent)
            sys.exit(1)
        
        # Combine things like cost, time, executable versions???
        parent_job['totalPrice'] = parent_job['totalPrice'] + child_job['totalPrice']
        parent_job['stoppedRunning'] += child_job.get('stoppedRunning') - child_job.get('startedRunning')
        parent_dxFile = dxencode.file_handler_from_fid(parent_fid)
        parent_job['step_parent_app_version'] = self.find_app_version(parent_dxFile)

        if verbose:
            print >> sys.stderr, "Step-child %s has step-parent %s" % (step_child,step_parent)
            print >> sys.stderr, json.dumps(parent_job,indent=4)
        return parent_job

    def make_payload_obj(self,out_type,rep_tech,fid,verbose=False):
        '''Returns an object for submitting a file to encode, with all dx info filled in.'''
        payload = {}
        payload['dataset'] = self.exp_id
        output_format = out_type.split('|')  # comes from hard-coded json key and *may* be "output_type|format|format_type"
        payload["output_type"] = output_format[0]  
        if len(output_format) > 2:
            payload["file_format_type"] = output_format[2]

        # Handle awkward condition where file object was posted, but then the file upload failed
        f_obj = self.found.get(fid)
        if f_obj != None:
            #assert f_obj['status'] == "upload failed"  # Actually this could be a posted file for which qc_metric post failed.
            payload['accession'] = f_obj['accession'] # accession tells dxencode.encoded_post_file() to patch obj, not post new
            payload['status'] = f_obj['status']

        dx_obj = dxencode.description_from_fid(fid)

        # Some steps are child processes which are not known by encodeD so we must proceed with the "step-parent"
        job = dxencode.job_from_fid(fid)
        if job.get('executableName') in self.STEP_CHILDREN.keys():
            child_job = job
            job = self.find_step_parent(fid,child_job,child_job['executableName'],self.STEP_CHILDREN[child_job['executableName']])

        if out_type != "QC_only": # Use new qc_object posting methods 
            payload["file_format"] = self.file_format(dx_obj["name"])
            if payload["file_format"] == None:
                print "Warning: file %s has unknown file format!" % dxencode.file_path_from_fid(fid)
            if len(output_format) > 1 and payload["file_format"] != output_format[1]:
                print "Warning: file %s has format %s but expecting %s!" % \
                                            (dxencode.file_path_from_fid(fid),payload["file_format"], output_format[1])
        payload["derived_from"] = self.find_derived_from(fid,job, verbose)
        payload['submitted_file_name'] = dxencode.file_path_from_fid(fid,projectToo=True)
        payload['file_size'] = dx_obj["size"]
        #payload['md5sum'] = calculated_md5 # Done i  validate_post applet
        if self.genome == None:
            print "ERROR: could not determine genome assembly! Check reference file properties."
            sys.exit(1)
        else:
            payload['assembly'] = self.genome
        if self.annotation != None:
            payload['genome_annotation'] = self.annotation
        elif self.exp_type in self.REQUIRE_ANNOTATION:
            print "ERROR: could not determine annotation on '%s' experiment! Check reference file properties." % (self.exp_type)
            sys.exit(1)
            

        dxFile = dxencode.file_handler_from_fid(fid)
        #versions = dxencode.get_sw_from_log(dxFile, '\* (\S+)\s+version:\s+(\S+)') # * STAR version: 2.4.0k
        versions = self.find_sw_versions(dxFile)
        notes = dxencode.create_notes(dxFile, versions)
        notes["notes_version"] = "5" # Cricket requests starting at "5", since earlier files uploads were distingusihed by user
        if 'totalPrice' in job:
            notes['dx_cost'] = "$" + str(round(job['totalPrice'],2))
        if 'step_parent_app_version' in job:
            notes['step_parent_app'] = job.get('step_parent_app_version') # One more little piece of the story.

        #print "  - Adding encoded information."
        payload = self.add_encoded_info(payload,rep_tech,fid)
        
        # Find or create step_run object
        step_run = self.enc_step_run_find_or_create(job,dxFile,rep_tech,test=self.test,verbose=verbose)
        #step_run = self.enc_step_run_find_or_create(job,dxFile,rep_tech,test=False,verbose=verbose)
        if step_run:
            if "dx_applet_details" in step_run:
                payload["step_run"] = "/analysis-step-runs/dnanexus:" + step_run["dx_applet_details"][0].get("dx_job_id")
            elif 'aliases' in step_run:
                payload["step_run"] = "/analysis-step-runs/" + dxencode.select_alias(step_run.get('aliases'))
            elif '@id' in step_run:
                payload["step_run"] = step_run.get('@id')
            #else: # What?
            if "notes" in step_run:
                step_run_notes = json.loads(step_run.get("notes"))
                # analysis_step and pipeline are calculated properties.
                #payload["analysis_step"] = "/analysis-steps/" + step_run_notes.get("step_name")
                #payload["pipeline"] = "/pipelines/encode:" + step_run_notes.get("pipeline_name")
                if "dx_analysis_id" in step_run_notes:
                    notes["dx_analysis_id"] = step_run_notes.get("dx_analysis_id")

        notes["dx_project_id"] = self.proj_id
        notes["dx_project_name"] = self.proj_name
        payload['notes'] = json.dumps(notes)
        
        if verbose:
            print >> sys.stderr, "payload:"
            print >> sys.stderr, json.dumps(payload,indent=4)
            print >> sys.stderr, "payload[notes]:"
            print >> sys.stderr, json.dumps(notes,indent=4)
        return payload


    def file_get_accession(self,fid,verify=False,fake_if_needed=False):
        '''Adds/replaces accession to a file's properties.'''
        # Posted files should have accession in properties
        fileDict = dxencode.description_from_fid(fid,properties=True)
        acc_key = dxencode.dx_property_accesion_key('https://www.encodeproject.org') # prefer production accession
        # check file properties
        accession = ''
        if "properties" in fileDict:
            if acc_key not in fileDict["properties"] and self.server_key != 'www':   # test, beta replicated from www
                acc_key = dxencode.dx_property_accesion_key(self.server)
            if acc_key in fileDict["properties"]:
                accession = fileDict["properties"][acc_key]
                # verify the file is posted
                if verify:
                    if fid not in self.found: 
                        f_obj =  self.enc_file_find_by_dxid(fid) # look by alias first:               
                        if f_obj == None:
                            f_obj = self.enc_lookup_json( 'files/' + accession,must_find=False) # look by accession
                        if f_obj != None: # Verifyably posted
                            self.found[fid] = f_obj
                        else:
                            accession = None
        
        if accession == None and fake_if_needed:
            if self.server_key == "www":
                accession = "ENCFF00FAKE"
            else:
                accession = "TSTFF00FAKE"
        return accession


    def file_mark_accession(self,fid,accession,test=True):
        '''Adds/replaces accession to a file's properties.'''
        acc_key = dxencode.dx_property_accesion_key(self.server)
        path = dxencode.file_path_from_fid(fid)
        acc = dxencode.dx_file_set_property(fid,acc_key,accession,add_only=True,test=test)
        if acc == None or acc != accession and not accession.endswith('FAKE'):
            print "Error: failed to update %s for file %s to '%s'" % (path,acc_key,accession)
        elif test:
            print "  - Test flag %s with %s='%s'" % (path,acc_key,accession)
        else:
            print "  - Flagged   %s with %s='%s'" % (path,acc_key,accession)


    def can_skip_validation(self,exp_type,payload,test=True):
        '''Returns True if validation can be skipped.'''
        if exp_type.startswith("long-rna-seq") or exp_type.startswith("small-rna-seq") or exp_type.startswith("rampage"):
            return True
        if payload["output_type"] in self.SKIP_VALIDATE.keys():
            return (self.SKIP_VALIDATE[payload["output_type"]] == payload["file_format"])
        return True

    def file_patch(self,fid,payload,test=True):
        '''Patches ENCODEd file with payload.'''
        f_obj = self.found.get(fid)
        if f_obj == None:
            return None
        accession = f_obj.get('accession')
        if accession == None:
            return None
            
        assert f_obj.get('accession') == payload.get('accession')
        assert f_obj.get('dataset') == "/experiments/" + payload.get('dataset') + "/"
        assert f_obj.get('file_format') == payload.get('file_format')
        assert f_obj.get('file_format_type') == payload.get('file_format_type')
        assert f_obj.get('output_type') == payload.get('output_type')
            
        update_payload = payload
        del update_payload['accession'] 
        del update_payload['status']
        del update_payload['dataset'] 
        del update_payload['file_format']
        del update_payload['file_format_type']
        del update_payload['output_type']

        if not test:
            ret = dxencode.encoded_patch_obj(accession, update_payload, self.server, self.authid, self.authpw)
            print "  * Patched '"+accession+"' with payload."
        else:
            print "  * Would patch '"+accession+"' with payload."
            #print json.dumps(update_payload,indent=4,sort_keys=True)
        return accession


    def file_post(self,fid,payload,test=True):
        '''Posts a file to encoded.'''
        path = payload['submitted_file_name'].split(':')[1]
        derived_count = len(payload["derived_from"])
        skip_validate = self.can_skip_validation(self.exp_type,payload)
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


    def print_analysis_totals(self,exp_id):
        '''Prints totals for all analyses for an experiment.'''
        if "ana_id" not in self.obj_cache["exp"]:
            return
        total_cost = 0
        total_dur = 0
        count = 0
        for ana_id in self.obj_cache["exp"]["ana_id"]:
            #try:
            ana = dxpy.api.analysis_describe(ana_id)
            total_cost += ana["totalPrice"]
            total_dur  += (ana["modified"] - ana["created"])
            count += 1
            #except:
            #    continue

        if total_cost > 0 or total_dur > 0:
            # 27888946 / 7.75 = 3598573.54838709677419
            duration = dxencode.format_duration(0,total_dur/1000,include_seconds=False)
            #   Print lrna.txt line as....  Then use grep ENC3 *.log | sed s/^.*\log://
            #print "cost:       GRCh38 v24 -     1,2   no     2016-02-17  2016-02-18  2016-02-18  2016-02-19  %s  $%.2f" % \
            #    (duration, total_cost)
            print "%s %d %s  cost: %s  $%.2f" % \
                (exp_id, len(self.obj_cache["exp"]["ana_id"]), self.obj_cache["exp"]["ana_id"][0], duration, total_cost)

    def run(self):
        '''Runs splasdown from start to finish using command line arguments.'''
        args = self.get_args()
        self.test = args.test
        self.qc_from_file = args.qc_from_file
        self.ignore = False
        self.genome = args.genome
        if args.ignore_properties:
            print "Ignoring DXFile properties (will post to test server)"
            self.ignore = args.ignore_properties
            self.server_key = 'test' # mandated because option is dangerous
        if args.way_back_machine:
            print "Using 'way back machine' to find files posted long ago."
            self.way_back_machine = args.way_back_machine
            
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
        total_patched = 0
        total_qc_objs = 0        
        for exp_id in args.experiments:
            sys.stdout.flush() # Slow running job should flush to piped log
            self.exp_id = exp_id
            self.obj_cache["exp"] = {}  # clear exp cache, which will hold exp specific wf_run and step_run objects
            # 1) Lookup experiment type from encoded, based on accession
            print "Working on %s..." % self.exp_id
            self.exp = dxencode.get_exp(self.exp_id,must_find=True,key=self.server_key)
            if self.exp == None or self.exp["status"] == "error":
                print "Unable to locate experiment %s in encoded (%s)" % (self.exp_id, self.server_key)
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
            files_to_post = self.find_needed_files(files_expected, test=self.test, verbose=args.verbose)
            print "- Found %d files that need to be handled" % len(files_to_post)
            if len(files_to_post) == 0:
                continue

            # 6) For each file that needs to be posted:
            exp_count += 1
            file_count = 0
            post_count = 0
            patch_count = 0
            qc_obj_count = 0
            for (out_type,rep_tech,fid,QC_only) in files_to_post:
                sys.stdout.flush() # Slow running job should flush to piped log
                file_name = dxencode.file_path_from_fid(fid)
                if args.start_at != None:
                    if not file_name.endswith(args.start_at):
                        continue
                    else:
                        print "- Starting at %s" % (file_name)
                        args.start_at = None
                    
                # a) discover all necessary dx information needed for post.
                # b) gather any other information necessary from dx and encoded.
                print "  Handle file %s" % dxencode.file_path_from_fid(fid)
                payload = self.make_payload_obj(out_type,rep_tech,fid, verbose=args.verbose)

                file_count += 1
                # c) Post file and update encoded database.
                if not QC_only: # Use new qc_object posting methods 
                    accession = self.file_patch(fid,payload,args.test) # rare cases a patch is all that is needed.
                    if accession != None:
                        patch_count += 1
                    else:
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
                        post_count += 1                        
                        self.file_mark_accession(fid,accession,args.test)  # This should have already been set by validate_post
                
                # e) After the file is posted, post or patch qc_metric objects.
                qc_count = self.handle_qc_metrics(fid,payload,test=self.test,verbose=args.verbose)
                if qc_count > 0:
                    print "  - handled %d qc_metric objects for file." % qc_count
                    qc_obj_count += 1
                    total_qc_objs += 1                                

                if args.files != 0 and file_count >= args.files:  # Short circuit for test
                    print "- Just trying %d file(s) by request" % file_count
                    break

            if not args.test:
                print "- For %s processed %d file(s), posted %d, patched %d, qc %d" % \
                                                        (self.exp_id, file_count, post_count, patch_count, qc_obj_count)
            else:
                print "- For %s processed %d file(s), would post %d, patched %d, qc %d" % \
                                                        (self.exp_id, file_count, post_count, patch_count, qc_obj_count)
            self.print_analysis_totals(self.exp_id)
            total_posted += post_count
            total_patched += patch_count

        if not args.test:
            print "Processed %d experiment(s), halted %d, posted %d file(s), patched %d file(s), %d qc object(s)" % \
                                                            (exp_count, halted, total_posted, total_patched, total_qc_objs)
        else:
            print "Processed %d experiment(s), halted %d, would post %d file(s), patched %d file(s), %d qc object(s)" % \
                                                            (exp_count, halted, total_posted, total_patched, total_qc_objs)
        if halted == exp_count:
            sys.exit(1)
        print "(finished)"


if __name__ == '__main__':
    '''Run from the command line.'''
    recovery = Splashdown()
    recovery.run()

