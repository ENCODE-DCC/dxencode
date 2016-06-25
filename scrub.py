#!/usr/bin/env python2.7
# scrub.py 1.0.0
#
# Scrub.py will remove all files for an experiment [replicate] and genome/annotation
#
# 1) Lookup experiment type from encoded, based on accession
# 2) Locate the experiment accession named folder
# 3) Given the experiment type, determine the expected results
# 4) Given expected results locate any files (by glob) that should be removed
#    a) each single replicate (in replicate sub-folders named as reN_N/
#    b) combined replicates in the experiment folder itself
# 5) For each file that should be removed, determine if the file has already been posted
# 6) For each file that needs to be removed and has already been posted, remove

import argparse,os, sys
import json, urlparse, subprocess, itertools, logging, time
from datetime import datetime
from base64 import b64encode
import commands

import dxpy
import dx
import encd

class Scrub(object):
    '''
    Scrub module removes posted experiment files from dnanexus.
    '''
    TOOL_IS = 'scrub'
    HELP_BANNER = "Scrubs posted files from DX. " 
    ''' This help banner is displayed by get_args.'''

    SERVER_DEFAULT = 'www'
    '''This the server to makes files have been posed to.'''

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
                                     "transcriptome alignments":                  "*_star_anno.bam",
                                     "QC_only":                                   "*_star_Log.final.out"         },
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
                "signals":         { "plus strand signal of all reads":           "*_srna_star_plusAll.bw",
                                     "minus strand signal of all reads":          "*_srna_star_minusAll.bw",
                                     "plus strand signal of unique reads":        "*_srna_star_plusUniq.bw",
                                     "minus strand signal of unique reads":       "*_srna_star_minusUniq.bw"     }  },
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
                "peaks":           { "transcription start sites|gff|gff3":        "*_rampage_peaks.gff.gz",
                                     "transcription start sites|bed|tss_peak":    "*_rampage_peaks.bed.gz",
                                     "transcription start sites|bigBed|tss_peak": "*_rampage_peaks.bb",
                                     "gene quantifications":                      "*_rampage_peaks_quant.tsv" } },
            "combined":   {
                "idr":             { "transcription start sites|bed|idr_peak":    "*_rampage_idr.bed.gz",
                                     "transcription start sites|bigBed|idr_peak": "*_rampage_idr.bb" },
                "mad_qc":          { "QC_only":                                   "*_mad_plot.png" }  },
        },
        "dna-me": {
            "step-order": [ "align","quantification","corr"], # How to: 1) combine 3 steps into 1; 2) tech lvl, bio lvl, exp lvl
            "replicate":  {
                "align":           { "alignments":  [ "*_techrep_bismark_pe.bam", "*_bismark.bam"  ] },   # *may* have samtools_flagstat, samtools_stats, Don't wan't bismark_map
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
    '''List of supported formats, and means of recognizing with file extensions.'''
    
    PRIMARY_INPUT_EXTENSION = [ "fastq","fq"]
    '''List of file extensions used to recognize primary inputs to parse accessions.'''
    
    def __init__(self):
        '''
        Scrub expects one or more experiment ids as arguments and will find files that 
        should be removed from the associated directory.
        '''
        self.args = {} # run time arguments
        self.server_key = 'www'  # TODO: replace with self.encd.server_key when Encd class is created
        self.server     = None    # TODO: replace with self.encd.server() when Encd class is created
        self.acc_prefix = "TSTFF"
        self.proj_name = None
        self.project = None
        self.proj_id = None
        self.exp = {}  # Will hold the encoded exp json
        self.exp_id = None
        self.exp_type = {}  # Will hold the experiment's assay_type, normalized to known tokens.
        self.genome = None  # genome should be required
        self.annotation = None  # if appropriate (mice), points the way to the sub-dir
        self.pipeline = None # pipeline definitions (filled in when experiment type is known)
        self.replicates = None # lost replicate folders currently found beneath experiment folder
        self.test = True # assume Test until told otherwise
        self.force = False # remove files whether posted or not
        self.remove_all = False # Removes experiment dir and all files beneath it recursively (Requires force!)
        logging.basicConfig(format='%(asctime)s  %(levelname)s: %(message)s')
        encd.logger = logging.getLogger(__name__ + '.dxe') # I need this to avoid some errors
        encd.logger.addHandler(logging.StreamHandler()) #logging.NullHandler)
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
                        help='One or more ENCODED experiment accessions, or file containing list',
                        nargs='+',
                        required=True)

        ap.add_argument('--project',
                        help="Project to run analysis in (default: '" + \
                                                         dx.env_get_current_project() + "')",
                        required=False)

        ap.add_argument('-f','--folder',
                        help="The location to search for experiment folders (default: " + \
                                                "'<project>:" + self.FOLDER_DEFAULT + "')",
                        default=self.FOLDER_DEFAULT,
                        required=False)

        ap.add_argument('--server',
                        help="Server that files should have been posted to (default: '" + self.SERVER_DEFAULT + "')",
                        default=self.SERVER_DEFAULT,
                        required=False)

        ap.add_argument('-g','--genome',
                        help="The genome assembly that files were aligned to (default: discovered if possible)",
                        default=None,
                        required=True)

        ap.add_argument('--test',
                        help='Test run only, do not launch anything.',
                        action='store_true',
                        required=False)

        ap.add_argument('--start_at',
                        help="Start processing with this file name (or possibly accession).",
                        default=None,
                        required=False)

        ap.add_argument('--files',
                        help="Just delete this number of files (default: all)",
                        type=int,
                        default=0,
                        required=False)

        ap.add_argument('--remove_all',
                        help='Remove all files and directory (default is to leave fastqs and workflows) Requires force!',
                        action='store_true',
                        required=False)

        ap.add_argument('--force',
                        help='Remove files regardless of whether they have been posted or not.',
                        action='store_true',
                        required=False)

        ap.add_argument('--verbose',
                        help='More debugging output.',
                        action='store_true',
                        required=False)

        if parse:
            return ap.parse_args()
        else:
            return ap


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


    def strip_comments(self,line,ws_too=False):
        """Strips comments from a line (and opptionally leading/trailing whitespace)."""
        bam = -1
        ix = 0
        while True:
            bam = line[ix:].find('#',bam + 1)
            if bam == -1:
                break
            bam = ix + bam
            if bam == 0:
                return ''
            if line[ bam - 1 ] != '\\':
                line = line[ 0:bam ]
                break  
            else: #if line[ bam - 1 ] == '\\': # ignore '#' and keep looking
                ix = bam + 1
                #line = line[ 0:bam - 1 ] + line[ bam: ]
                
        if ws_too:
            line = line.strip()
        return line 


    def load_exp_list(self,exp_ids,verbose=False):
        '''Returns a sorted list of experiment accessions from command-line args.'''
        #verbose=True
        id_list = []
        file_of_ids = None
        
        # If only one, it could be a file
        if len(exp_ids) == 1:
            candidate = exp_ids[0]
            if candidate.startswith("ENCSR") and len(candidate) == 11:
                id_list.append(candidate)
                return id_list
            else:
                file_of_ids = candidate
                if file_of_ids != None:
                    with open(file_of_ids, 'r') as fh:
                        #line = fh.readline()
                        for line in fh:
                            line = self.strip_comments(line,True)
                            if line == '':
                                continue
                            candidate = line.split()[0]
                            if candidate.startswith("ENCSR") and len(candidate) == 11:
                                id_list.append(candidate)
                            elif verbose:
                                print >> sys.stderr, "Value is not experiment id: '"+candidate+"'"
                    
        elif len(exp_ids) > 0:
            for candidate in exp_ids:
                if candidate.startswith("ENCSR") and len(candidate) == 11:
                    id_list.append(candidate)
                elif verbose:
                    print >> sys.stderr, "Value is not experiment id: '"+candidate+"'"
        
        if len(id_list) > 0:
            sorted_exp_ids = sorted(id_list)
            if verbose:
                print >> sys.stderr, "Experiment ids: "
                print >> sys.stderr, json.dumps(sorted_exp_ids)
            print "Requested scrubbing %d experiments" % len(sorted_exp_ids)
            return sorted_exp_ids
            
        return []


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


    def find_step_files(self,file_globs,result_folder,rep_tech,verbose=False):
        '''Returns tuple list of (type,rep_tech,fid) of ALL files expected for a single step.'''
        step_files = []

        for token in file_globs.keys():
            if type(file_globs[token]) == list:
                for file_glob in file_globs[token]:
                    if token != "QC_only": 
                        if self.file_format(file_glob) == None:
                            print "Error: file glob %s has unknown file format! Please fix" % file_glob
                            sys.exit(1)
                    if verbose:
                        print >> sys.stderr, "-- Looking for %s" % (result_folder + file_glob)
                    fid = dx.find_file(result_folder + file_glob,self.proj_id, recurse=False)
                    if fid != None:
                        QC_only = (token == "QC_only") # Use new qc_object posting methods  
                        step_files.append( (token,rep_tech,fid,QC_only) )
                        break # Only looking for the first hit               
            else:
                if token != "QC_only": 
                    if self.file_format(file_globs[token]) == None:
                        print "Error: file glob %s has unknown file format! Please fix" % file_globs[token]
                        sys.exit(1)
                if verbose:
                    print >> sys.stderr, "-- Looking for %s" % (result_folder + file_globs[token])
                fid = dx.find_file(result_folder + file_globs[token],self.proj_id, recurse=False)
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


    def input_exception(self,inp_fid):
        '''Returns True if this is one of a limit number of input files we do not track in encodeD.'''
        # TODO: move specifics to json at top of file.
        # Unfortunate special case: the map_report is essentially a QC_only file but is an input to a step in order to 
        # combine multiple map_reports into a single qc_metric.
        try:
            if self.exp_type != "dna-me" or not dx.file_path_from_fid(inp_fid).endswith("_map_report.txt"):
                #print "** Ignoring file: " + dx.file_path_from_fid(inp_fid)
                return False
        except:
            pass
        return True
                
    def find_removable_files(self,files_expected,test=True,verbose=False):
        '''Returns the tuple list of files that NEED to be posted to ENCODE.'''
        removable = []
        not_posted = 0
        acc_key = dx.property_accesion_key(self.server_key) # 'accession'

        for (out_type, rep_tech, fid, QC_only) in files_expected:
            if not QC_only:
                acc = dx.file_get_property(acc_key,fid)
                if acc != None or self.force:
                    removable.append( (out_type,rep_tech,fid, False) )
                elif self.input_exception(fid):
                    removable.append( (out_type,rep_tech,fid, False) )
                else:
                    # TODO: back up plan, look on encodeD?
                    not_posted += 1
                    print >> sys.stderr, "* WARNING: file '" + dx.file_path_from_fid(fid) + \
                                                                             "' has not been posted, and will not be scrubbed."
            #else:  # TODO: How to handle qc_only files (other than just removing them)?

        if not_posted > 0 and not self.force: # If even one file is not posted, then none are removable
            return []
            
        # if all expected non-QC files are remobable, then go ahead and remove the qc ones as well
        for (out_type, rep_tech, fid, QC_only) in files_expected:
            if QC_only:
                removable.append( (out_type,rep_tech,fid, True) )
        
        if verbose:
            print >> sys.stderr, "Removable files:"
            print >> sys.stderr, json.dumps(removable,indent=4)
        return removable        

    def run(self):
        '''Runs scrub from start to finish using command line arguments.'''
        args = self.get_args()
        self.test = args.test
        self.genome = args.genome
        self.force = args.force
        self.remove_all = args.remove_all
            
        self.server_key = args.server
        encd.set_server_key(self.server_key) # TODO: change to self.encd = Encd(self.server_key)
        self.server = encd.get_server()        
        
        if self.server_key == "www":
            self.acc_prefix = "ENCFF"
        self.proj_name = dx.env_get_current_project()
        if self.proj_name == None or args.project != None:
            self.proj_name = args.project
        if self.proj_name == None:
            print "Please enter a '--project' to run in."
            sys.exit(1)

        self.project = dx.get_project(self.proj_name)
        self.proj_id = self.project.get_id()
        print "== Running in project [%s] and expect files already posted to the [%s] server ==" % \
                                                        (self.proj_name,self.server_key)

        self.exp_ids = self.load_exp_list(args.experiments,verbose=args.verbose)
        if len(self.exp_ids) == 0:
            print >> sys.stderr, "No experiment id's requested."
            self.ap.print_help()
            sys.exit(1)       

        exp_count = 0
        exp_removed = 0
        exp_kept = 0
        deprecates_removed = 0
        total_removed = 0
        for exp_id in self.exp_ids:
            dx.clear_cache()
            sys.stdout.flush() # Slow running job should flush to piped log
            self.exp_id = exp_id
            # 1) Lookup experiment type from encoded, based on accession
            print "Working on %s..." % self.exp_id
            self.exp = encd.get_exp(self.exp_id,must_find=True)
            if self.exp == None or self.exp["status"] == "error":
                print "Unable to locate experiment %s in encoded (%s)" % (self.exp_id, self.server_key)
                continue
            self.exp_type = encd.get_exp_type(self.exp_id,self.exp,self.EXPERIMENT_TYPES_SUPPORTED)
            if self.exp_type == None:
                continue

            # 2) Locate the experiment accession named folder
            # NOTE: genome and annotation are not known for this exp yet, so the umbrella folder is just based on exp_type
            self.umbrella_folder = dx.umbrella_folder(args.folder,self.FOLDER_DEFAULT,self.proj_name,self.exp_type,"posted",self.genome)
            if args.test:
                print "- Umbrella folder: " + self.umbrella_folder
                
            self.exp_folder = dx.find_exp_folder(self.project,exp_id,self.umbrella_folder,warn=True)
            if self.exp_folder == None:
                continue
            exp_count += 1
            print "- Examining %s:%s for '%s' results..." % \
                                            (self.proj_name, self.exp_folder, self.exp_type)

            # Could be quick... remove everything!
            if self.remove_all and self.force:
                exp_removed += 1
                if self.test: 
                    print "* Would remove %s:%s and all results within..." % (self.proj_name, self.exp_folder)
                else:
                    print "* Removing %s:%s and all results within..." % (self.proj_name, self.exp_folder)
                    dxpy.api.project_remove_folder(self.proj_id,{'folder':self.exp_folder,'recurse':True})
                continue
            
            # Remove any 'deprecated' subfolder
            deprecated_folder = self.exp_folder + "deprecated/"
            if dx.project_has_folder(self.project, deprecated_folder): 
                deprecates_removed += 1
                if self.test: 
                    print "* Would remove %s:%s and all results within..." % (self.proj_name, deprecated_folder)
                else:
                    print "* Removing %s:%s and all results within..." % (self.proj_name, deprecated_folder)
                    dxpy.api.project_remove_folder(self.proj_id,{'folder':deprecated_folder,'recurse':True})

                       
            # 3) Given the experiment type, determine the expected results
            self.pipeline   = self.pipeline_specification(args,self.exp_type,self.exp_folder)
            self.replicates = dx.find_replicate_folders(self.project,self.exp_folder, verbose=args.verbose)

            # 4) Given expected results locate any files (by glob) that should have been posted for
            #    a) each single replicate (in replicate sub-folders named as reN_N/
            #    b) combined replicates in the experiment folder itself
            files_expected = self.find_expected_files(self.exp_folder, self.replicates, verbose=args.verbose)
            print "- Found %d files that are available to remove." % len(files_expected)
            if len(files_expected) == 0:
                continue

            # 5) For each file that is available to be removed, determine if the file has been posted first.
            files_to_remove = self.find_removable_files(files_expected, test=self.test, verbose=args.verbose)
            print "- Found %d files that may be removed" % len(files_to_remove)
            if len(files_to_remove) == 0:
                print "- KEEPING: If even one file has not been posted, no files may be removed without force."
                exp_kept += 1
                continue

            # 6) For each file that needs to be removed:
            files_removed = 0
            for (out_type,rep_tech,fid,QC_only) in files_to_remove:
                sys.stdout.flush() # Slow running job should flush to piped log

                if args.files != 0 and file_count >= args.files:  # Short circuit for test
                    print "- Just trying %d file(s) by request" % file_count
                    partial = True
                    break

                file_name = dx.file_path_from_fid(fid)
                if args.start_at != None:
                    if not file_name.endswith(args.start_at):
                        continue
                    else:
                        print "- Starting at %s" % (file_name)
                        args.start_at = None
                    
                if self.test:
                    print "  * Would remove file %s..." % file_name
                else:
                    print "  * Removing file %s..." % file_name
                    dxpy.api.project_remove_objects(self.proj_id,{'objects':[fid]})
                files_removed += 1
                
            if not args.test:
                print "- For %s processed %d file(s), removed %d files" % (self.exp_id, len(files_expected), files_removed)
            else:
                print "- For %s processed %d file(s), would remove %d files" % (self.exp_id, len(files_expected), files_removed)
            total_removed += files_removed
            
        if not args.test:
            print "Processed %d experiment(s), erased %d, kept %d, removed %d deprecate folder(s) and %d file(s)" % \
                                                      (exp_count, exp_removed, exp_kept, deprecates_removed, total_removed)
        else:
            print "Processed %d experiment(s), would erase %d, would keep %d, would remove %d deprecate folder(s) and %d file(s)" % \
                                                      (exp_count, exp_removed, exp_kept, deprecates_removed, total_removed)
        print "(finished)"


if __name__ == '__main__':
    '''Run from the command line.'''
    scrub = Scrub()
    scrub.run()

