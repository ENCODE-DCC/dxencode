#!/usr/bin/env python2.7
# mission_log.py 0.1.0
#
# Mission_log is meant to report on experiments that have been posted to encoded.
# It expects a list of experiment accessions to report on.
#
# Basic design:
# For each experiment:
# a) find relevant files in determined order 
# b) find metric associated (quality_metric object or other):
# c) print 1 line per replicate and requested fields per metric. 
# d) print 1 combined line for any combined metrics
# e) rinse and repeat
# Print totals
#
#star: "Number of input reads"
#      "Uniquely mapped reads number" 
##exp        rep1_reads rep1_mapped_unique rep1_u_pct rep2_reads rep2 mapped_unique rep2_u_pct        mad pearson spearman
#ENCSR717SJA 23943208   20967172           87.57%     24354505   21293182           87.43%     0.176  0.9883518   0.9911513

#end of report
#max
#min
#avg

import argparse,os, sys
import json, urlparse, subprocess, itertools, logging, time
#import requests, re, shlex, time
from datetime import datetime
from base64 import b64encode
import commands, string

import dxpy
import dx
import encd

class Mission_log(object):
    '''
    Mission_log module reports information for a set of experiments posted to ENCODEd.
    '''
    TOOL_IS = 'mission_log'
    HELP_BANNER = "Reports statistics for a set of experiments posted to ENCODEd. " 
    ''' This help banner is displayed by get_args.'''

    SERVER_DEFAULT = 'www'
    '''This the default server to report from.'''

    EXPERIMENT_TYPES_SUPPORTED = [ 'long-rna-seq', 'small-rna-seq', 'rampage', 'dna-me'] #, 'rampage','dnase','dna-me','chip-seq' ]
    '''This module supports only these experiment (pipeline) types.'''

    REPORT_DEFAULT = 'star-mad'
    '''This the default report type.'''

    ASSEMBLIES_SUPPORTED = { "hg19": "hg19", "hg38": "GRCh37", "mm10": "mm10" }
    '''This module supports only these assemblies.'''

    ANNOTATIONS_SUPPORTED = [ 'V24', 'V19', 'M4' ]
    '''This module supports only these annotations.'''
    
    FOLDER_DEFAULT = "/"
    '''Where to start the search for experiment folders in DX.'''

    # REPORT_SPECS is a dict covering all supported types of reports.
    # <report_type>: A dict keyed on report_type containing report specific json.
    #   "sources": a list of one or both of 'encodeD' or 'DX'  Not all reports can use both sources for content.
    #   <pipeline or assay>: A dictionary keyed on the assay_type that furter describes what to look for.
    #      "output_types": list of encodeD FILE output_types that are needed to complete the report.
    #      <out_type>: Dict keyed on either output_type or sub_type further defining the files and metrics to be gathered
    #         "sub_types": list of sub_types for an output_type if they exist.
    #         "suffix":    List of one or more file name endings, used to find files in DX or verify files in encodeD
    #         "metrics":   List of one or more metric keys to the METRICS DEF defined below.   
    REPORT_SPECS = {
         'star-mad': {
             "sources":       ["encodeD","DX"],
             "long-rna-seq":  { "output_types": [ "transcriptome alignments", "gene quantifications" ],
                                "transcriptome alignments": { "suffix":  [ "_star_anno.bam" ],
                                                              "metrics": [ "star" ] }, 
                                "gene quantifications":     { "suffix":  [ "_rsem.genes.results", "_mad_plot.png" ],
                                                              "metrics": [ "mad" ] },
                              },
             "small-rna-seq": { "output_types": [ "alignments", "gene quantifications" ],
                                "alignments":           { "suffix":  [ "_srna_star.bam" ],
                                                          "metrics": [ "star" ] }, 
                                "gene quantifications": { "suffix":  [ "_srna_star_quant.tsv", "_mad_plot.png" ],
                                                          "metrics": [ "mad" ] },
                              },
             "rampage":       { "output_types": [ "alignments", "gene quantifications" ],
                                "alignments":           { "suffix":  [ "_rampage_star_marked.bam" ],
                                                          "metrics": [ "star" ] }, 
                                "gene quantifications": { "suffix":  [ "_rampage_peaks_quant.tsv", "_mad_plot.png" ],
                                                          "metrics": [ "mad" ] },
                              }
         },
         'tophat_times': {
             "sources":       ["encodeD","DX"],
             "long-rna-seq":  { "output_types": [ "alignments" ],
                                "alignments":   { "subtypes": [ "star", "tophat" ] },
                                "star":         { "suffix":   [ "_star_genome.bam" ],
                                                  "metrics":  [ "star2", "file_cost" ] }, 
                                "tophat":       { "suffix":   [ "_tophat.bam" ], 
                                                  "metrics":  [ "file_cost" ] }, 
                              }
         },
         'wgbs': {
             "sources":       ["encodeD","DX"],
             "dna-me":        { #"output_types": [ "alignments","methylation state at CpG" ],
                                #"alignments":   { "suffix":   [ "_techrep_bismark_pe.bam", "_techrep_bismark_se.bam" ],
                                #                  "metrics":  [ "bismark", "samtools_flagstats" ] },
                                "output_types": [ "methylation state at CpG" ], 
                                "methylation state at CpG": { 
                                                  "suffix":   [ "_biorep_CpG.bed.gz", "_CpG_corr.txt" ], 
                                                  "metrics":  [ "samtools_flagstats", "bismark", "cpg_correlation" ] }, 
                              }
         },
         'cost': {
             "sources":       ["encodeD","DX"], # Unfortunately this is only of DX.
             "long-rna-seq":  { "output_types": [ "gene quantifications", "mad-qc" ],
                                "gene quantifications": { "suffix":  [ "_rsem.genes.results" ],
                                                          "metrics": [ "rep_cost" ] },
                                "mad-qc": { "suffix":  [ "_mad_plot.png" ],
                                            "metrics": [ "exp_cost" ] }, 
                              },
             "small-rna-seq": { "output_types": [ "signal", "mad-qc" ],
                                "signal": { "suffix":  [ "_minusUniq.bw" ],
                                            "metrics": [ "rep_cost" ] },
                                "mad-qc": { "suffix":  [ "_mad_plot.png" ],
                                            "metrics": [ "exp_cost" ] }, 
                              },
             "rampage":       { "output_types": [ "gene quantifications", "idr" ],
                                "gene quantifications": { "suffix":  [ "_rampage_peaks_quant.tsv", "_mad_plot.png" ],
                                                          "metrics": [ "rep_cost" ] },
                                "idr": { "suffix":  [ "_rampage_idr.png" ],
                                            "metrics": [ "exp_cost" ] }, 
                              },
             "dna-me":        { "output_types": [ "methylation state at CpG", "correlation" ], 
                                "methylation state at CpG": { 
                                                  "suffix":   [ "_biorep_CpG.bed.gz" ], 
                                                  "metrics":  [ "rep_cost" ] }, 
                                "correlation": { 
                                                  "suffix":   [ "_CpG_corr.txt" ], 
                                                  "metrics":  [ "exp_cost" ] }, 
                              }
         },
    }
    '''For each report type, these are the mappings for file 'output_types' to 'quality_metric' types.'''
    
    # METRIC_DEFS is a dict describing the metrics that are collected, one or several are needed for each exp/file in a report.
    # <metric_type>: A dict describing a single metric object.
    #   "per": contents of a metric object pertain to either a 'replicate' or the whole 'experiment'
    #   "name": (optional, default: metric_type) used to find the qc_metric in encodeD 
    #                             (e.g samtools_flagstats will be used to find 'SamtoolsFlagstatsQualityMetric')
    #   "dx_key": (optional, default: metric_type) Used to find the metric in a DX file's details json blob.
    #   "columns": An ordered  list of the properties within a quality_metric object to include in the report.
    #   "headings": (optional, default: columns). A mapping of column properties to column headings.
    #   "special_metric": (optional, default: False). Create metric from non-quality_metric object info.
    #   "prepend": (optional) If true, prepend the out_type from REPORT_SPECS to the column headers (e.g. "Star Cost").
    #   "format": (optional) Special formatting instructions.
    METRIC_DEFS = {
       'star': {"per": "replicate",
                "dx_key": "STAR_log_final",
                "columns": [
                            "Number of input reads",
                            "Uniquely mapped reads number",
                            "Uniquely mapped reads %" ],
                "headings": {
                                "Number of input reads":        "Reads",
                                "Uniquely mapped reads number": "Unique Reads",
                                "Uniquely mapped reads %":      "Unique Read %"}
              },
       'star2': {"per": "replicate",
                "name": "star",
                "dx_key": "STAR_log_final",
                "columns": [
                            "Number of input reads",
                            "Uniquely mapped reads %",
                            "% of reads unmapped: too short" ],
                "headings": {
                                "Number of input reads":         "Reads",
                                "Uniquely mapped reads %":       "Unique %",
                                "% of reads unmapped: too short":"Unmapped %"}
              },
       'mad': { "per": "experiment",
                "dx_key": "MAD.R",
                "columns": [
                            "MAD of log ratios",
                            "Pearson correlation",
                            "Spearman correlation" ],
              },
       'file_cost': { "per": "replicate",
                      "special_metric": True,
                      "prepend": "type",
                      "columns": [
                                    "File Size",
                                    "Time",
                                    "Cost" ],
                      "format": { "File Size": "file_size",
                                  "Time": "duration" },
                      },
         # TODO: Make a total time/cost metric
       'rep_cost': {  "per": "replicate",
                      "special_metric": True,
                      "columns": [ "CPU Time","Job Time","Total Cost" ],
                      "format": { "CPU Time": "duration", "Job Time": "duration" },
                      },
       'exp_cost': {  "per": "experiment",
                      "special_metric": True,
                      "columns": [ "CPU Time","Job Time","Total Cost" ],
                      "sumable": [ "CPU Time","Job Time","Total Cost" ],
                      "format": {  "CPU Time": "duration", "Job Time": "duration" },
                      },
       'samtools_flagstats': { 
                "per": "replicate",
                "columns": [
                            "mapped",
                            "mapped_pct" ],
                "headings": {
                            "mapped":         "Mapped",
                            "mapped_pct":     "Mapped %"}
              },
       'bismark': { 
                "per": "replicate",
                "dx_key": "bismark_map",
                "columns": [
                            "C methylated in CpG context",
                            "C methylated in CHG context",
                            "C methylated in CHH context",
                            "lambda C methylated in CpG context",
                            "lambda C methylated in CHG context",
                            "lambda C methylated in CHH context",
                            "Mapping efficiency" ],
                "headings": {
                            "C methylated in CpG context":         "Methylated: CpG",
                            "C methylated in CHG context":         "CHG",
                            "C methylated in CHH context":         "CHH",
                            "lambda C methylated in CpG context":  "Lambda Methylated: CpG",
                            "lambda C methylated in CHG context":  "CHG",
                            "lambda C methylated in CHH context":  "CHH",
                            "Mapping efficiency":                  "Efficiency" }
              },
       'cpg_correlation': { 
                "per": "experiment",
                "dx_key": "bedmethyl_corr",
                "columns": [
                            "CpG pairs",
                            "CpG pairs with atleast 10 reads each",
                            "Pearson Correlation Coefficient" ],
                "headings": {
                            "CpG pairs with atleast 10 reads each":  "with >= 10 reads",
                            "Pearson Correlation Coefficient":      "Pearson" }
              },
    }

    def __init__(self):
        '''
        Mission_log generates reports from pipeline results that include quality_metrics that have already been 
        uploaded to encodeD.  
        '''
        self.args = {} # run time arguments
        self.server_key = 'www'  # TODO: replace with self.encd.server_key when Encd class is created
        self.exp = {}  # Will hold the encoded exp json
        self.exp_id = None
        self.exp_type = None  # Will hold the experiment's assay_type, normalized to known tokens.
        self.genome = None
        self.genome_warning = False
        self.report_specs = {} # Unless the proper report specifications are discovered, a report can't be run
        self.data_mine       = None # Where to get QC metrics and other data from: 'encodeD' or 'DX'
        self.proj_name       = None # Only needed when qc_source is 'DX'
        self.project         = None # Only needed when qc_source is 'DX'
        self.umbrella_folder = None # Only needed when qc_source is 'DX'
        self.folder          = None # Where the user says to start looking for files
        self.obj_cache = {} # certain things take time to find or create and are needed multiple times
        print >> sys.stderr, " "


    def get_args(self,parse=True):
        '''Parse the input arguments.'''
        ### PIPELINE SPECIFIC
        ap = argparse.ArgumentParser(description=self.HELP_BANNER + "All report texts " +
                    "is written to stdout, while progress and any warnings or errors is written to stderr. ")
        ### PIPELINE SPECIFIC

        ap.add_argument('-e', '--experiments',
                        help='One or more ENCODED experiment accessions',
                        nargs='+',
                        required=False)

        ap.add_argument('-f','--file',
                        help="A file containing experiment accession (first word in line, # comments ignored)",
                        default=None,
                        required=False)

        ap.add_argument('-r','--report-type',
                        help="The report type to print (supported: "+str(self.REPORT_SPECS.keys())+") " + \
                                                      "(default: '" + self.REPORT_DEFAULT + "')",
                        default=self.REPORT_DEFAULT,
                        required=False)

        ap.add_argument('-g','--genome',
                        help="The genome assembly to run on (default: '" + \
                                                                  dx.GENOME_DEFAULTS['human'] + "')",
                        default=None,
                        required=False)

        ap.add_argument('-d','--dx',
                        help="Get content from 'DX' instead of the default 'encodeD'.",
                        action='store_true',
                        required=False)

        ap.add_argument('--folder',
                        help="If getting content from DX, the location to search for experiment folders (default: " + \
                                                "'<project>:" + self.FOLDER_DEFAULT + "')",
                        default=self.FOLDER_DEFAULT,
                        required=False)

        ap.add_argument('--server',
                        help="Server on which to look for posted files (default: '" + self.SERVER_DEFAULT + "')",
                        default=self.SERVER_DEFAULT,
                        required=False)

        ap.add_argument('--verbose',
                        help='More debugging output.',
                        action='store_true',
                        required=False)

        if parse:
            self.ap = ap
            return ap.parse_args()
        else:
            return ap


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


    def load_exp_list(self,exp_ids,file_of_ids,verbose=False):
        '''Returns a sorted list of experiment accessions from command-line args.'''
        #verbose=True
        id_list = []
        if exp_ids != None and len(exp_ids) > 0:
            for candidate in exp_ids:
                if candidate.startswith("ENCSR") and len(candidate) == 11:
                    id_list.append(candidate)
                elif verbose:
                    print >> sys.stderr, "Value is not experiment id: '"+candidate+"'"
        
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
                    
        if len(id_list) > 0:
            sorted_exp_ids = sorted(id_list)
            if verbose:
                print >> sys.stderr, "Experiment ids: "
                print >> sys.stderr, json.dumps(sorted_exp_ids)
            return sorted_exp_ids
            
        return []


    def find_replicates(self, exp_id, exp, verbose=False):
        '''Returns a list of replicates with input files for this experiment.'''
            
        # Must look through exp and find all replicates!
        if exp != self.exp or "full_mapping" not in self.obj_cache["exp"]:
            self.obj_cache["exp"]["full_mapping"] = encd.get_full_mapping(exp_id,exp)
        replicates = encd.get_reps(exp_id, exp=exp, full_mapping=self.obj_cache["exp"]["full_mapping"])
        for rep in replicates:
            if self.genome == None:
                if rep['organism'] in dx.GENOME_DEFAULTS:
                    self.genome = dx.GENOME_DEFAULTS[rep['organism']]
                else:
                    print >> sys.stderr, "Organism %s not currently supported" % rep['organism']
                    sys.exit(1)
            elif self.genome != dx.GENOME_DEFAULTS[rep['organism']]:
                if self.genome_warning:
                    print >> sys.stderr, "WARNING: Mixing genomes in one report: %s and %s.  Will only report '%s'" % \
                                                        (self.genome, dx.GENOME_DEFAULTS[rep['organism']],self.genome)
                    self.genome_warning  = False
                
        if verbose:
            print >> sys.stderr, "Replicates:"
            print >> sys.stderr, json.dumps(replicates,indent=4)
        return replicates


    def get_relevant_dx_files(self, exp_id, exp, report_specs, verbose=False):
        '''Returns list of DX files for selected files available.'''
        #verbose=True
        dx_files = []
        self.obj_cache["exp"]["files"] = {}
        
        # Need a project specific place to start looking for files
        if self.umbrella_folder == None:
            self.umbrella_folder = dx.umbrella_folder(self.folder,self.FOLDER_DEFAULT,self.proj_name,self.exp_type)
            if self.folder == self.FOLDER_DEFAULT:
                self.umbrella_folder = self.umbrella_folder.replace("/runs/","/posted/")
            if self.genome not in self.umbrella_folder:
                self.umbrella_folder = self.umbrella_folder + self.genome + '/'
            if verbose:
                print >> sys.stderr, "Using umbrella folder '%s'" % (self.umbrella_folder)
            
        # Find experiment dir
        if verbose:
            print >> sys.stderr, "Looking in '%s:%s' for folder for %s." % (self.proj_name, self.umbrella_folder, exp_id)
        exp_folder = dx.find_exp_folder(self.project,exp_id,self.umbrella_folder)
        if exp_folder == None:
            print >> sys.stderr, "ERROR: Can't find experiment folder for %s underneath '%s:%s.'" % \
                                        (exp_id, self.proj_name, self.umbrella_folder)
            return dx_files
        elif verbose:
            print >> sys.stderr, "- Examining %s:%s for '%s' results..." % \
                                        (self.proj_name, exp_folder, self.exp_type)
        # Now look for replicate dirs:
        rep_folders = dx.find_replicate_folders(self.project,exp_folder)
        if len(rep_folders) == 0:
            print >> sys.stderr, "ERROR: Can't find any replicate folders in %s." % (self.proj_name, exp_folder)
            return dx_files
        elif verbose:
            print >> sys.stderr, "- Found %d replicate folders." % len(rep_folders)

        # Now we can look for DX files:
        fids_found = []
        for rep_folder in rep_folders:
            if verbose:
                print >> sys.stderr, "- Looking in folder %s." % rep_folder
            for out_type in report_specs["output_types"]:
                types_list = [ out_type ] 
                if "subtypes" in report_specs[out_type]:
                    types_list = report_specs[out_type]["subtypes"] 
                for sub_type in types_list:
                    for suffix in report_specs[sub_type]["suffix"]:
                        if verbose:
                            print >> sys.stderr, "- Looking for files for type '%s' with suffix '%s'." % (sub_type,suffix)
                        fid = dx.find_file(exp_folder + rep_folder + '/*' + suffix,self.proj_name, recurse=False)
                        if fid == None: # Or check in experiment folder
                            if verbose:
                                print >> sys.stderr, "- Looking in experiment folder %s for files with suffix '%s'." % \
                                                                                                            (exp_folder,suffix)
                            fid = dx.find_file(exp_folder + '*' + suffix,self.proj_name, recurse=False)
                        if fid == None:
                            continue
                        if fid in fids_found: # Did we already get this file at the experiment level?
                            continue
                        fids_found.append(fid)
                        if verbose:
                            print >> sys.stderr, "- Found file '%s'." % fid
                        file_dx_obj = dx.description_from_fid(fid,properties=True)
                        qc_json = dx.file_get_details(fid)
                        if qc_json and "QC" in qc_json:  # Not likely but QC json could be subsection of details
                            qc_json = qc_json["QC"]
                        if not qc_json:
                            # If this report type uses "special metrics" then there may not be any DX qc metric.
                            special = False
                            for metric_id in report_specs[sub_type]["metrics"]:
                                if self.METRIC_DEFS[metric_id].get("special_metric",False):
                                    special = True
                                    break
                            if not special:
                                if verbose:
                                    print >> sys.stderr, "- Found file '%s' has no qc_josn, so skipping it." % \
                                                                                                (file_dx_obj['name'])
                                continue
                        else:
                            if verbose:
                                print >> sys.stderr, "- Found qc blob for %s:" % fid
                                #print >> sys.stderr, json.dumps(qc_json,indent=4,sort_keys=True)
                            file_dx_obj["QC"] = qc_json    
                        self.obj_cache["exp"]["files"][fid] = file_dx_obj
                        dx_files.append((out_type,sub_type,rep_folder,fid))
        
        if verbose:
            print >> sys.stderr, "DX files: %d" % len(dx_files)
            for (out_type,sub_type,rep_folder,fid) in dx_files:
                print >> sys.stderr, "'%s' %s %s %s" % \
                        (out_type,rep_folder,fid,self.obj_cache["exp"]["files"][fid]['name'])
        return dx_files

 
    def find_matching_enc_files(self, files, out_type, rep_tech, sub_type=None, suffix=None):
        '''Returns list of matched files for out_type, replicate and possible ending.'''
        
        matching_files = []
        for file_enc_obj in files:
            if file_enc_obj.get('output_type') != out_type:
                continue
            if self.genome != None and file_enc_obj.get('assembly',self.genome) != self.genome:
                continue
            if suffix != None:
                file_name = file_enc_obj["submitted_file_name"]
                if not file_name.endswith(suffix):
                    continue
            if 'replicate' not in file_enc_obj:
                print >> sys.stderr, "WARNING: %s %s has no 'replicate'" % \
                                    (file_enc_obj['accession'],file_enc_obj['submitted_file_name']) 
                continue
            br = file_enc_obj['replicate']['biological_replicate_number']
            tr = file_enc_obj['replicate']['technical_replicate_number']
            file_rep_tech = "rep%d_%d" % (br,tr)
            if file_rep_tech != rep_tech:
                continue
            acc = file_enc_obj['accession'] 
            matching_files.append((out_type,sub_type,rep_tech,acc))
            self.obj_cache["exp"]["files"][acc] = file_enc_obj
        return matching_files
        

    def get_relevant_enc_files(self, exp_id, exp, report_specs, reps, verbose=False):
        '''Returns list of encodeD file objects for selected files available on encoded.'''
        #verbose=True

        enc_files = []
        self.obj_cache["exp"]["files"] = {}
        
        files = encd.get_exp_files(exp,report_specs["output_types"],lab="encode-processing-pipeline")
        # special case to get around m2,m3 files
        new_list = []
        while len(files) > 0:
            file_enc_obj = files.pop()
            if self.genome != None and file_enc_obj.get('assembly',self.genome) != self.genome:
                continue
            if "genome_annotation" not in file_enc_obj or file_enc_obj["genome_annotation"] in self.ANNOTATIONS_SUPPORTED:
               new_list.append(file_enc_obj)
        files = new_list
        if len(files) == 0:
            if verbose:
                print >> sys.stderr, "Experiment "+exp_id+" has no results to report."
            return files

        # How to place them in order?
        for rep in reps:
            for out_type in report_specs["output_types"]:
                types_list = [ out_type ] 
                if "subtypes" in report_specs[out_type]:
                    types_list = report_specs[out_type]["subtypes"] 
                for sub_type in types_list:
                    for suffix in report_specs[sub_type]["suffix"]:
                        matching_files = self.find_matching_enc_files(files, out_type, rep['rep_tech'], sub_type, suffix)
                        if len(matching_files) > 0:
                            enc_files.extend(matching_files)
    
        if verbose:
            print >> sys.stderr, "Encoded files:"
            for (out_type,sub_type, rep_tech,acc) in enc_files:
                print >> sys.stderr, "'%s' %s %s %s" % \
                        (out_type,rep_tech,acc,self.obj_cache["exp"]["files"][acc]['submitted_file_name'])
        return enc_files


    def get_relevant_files(self, exp_id, exp, report_specs, reps, verbose=False):
        '''Returns list of encodeD or DX file objects for selected files available on encoded.'''

        if self.data_mine == "DX":
            return self.get_relevant_dx_files(exp_id, exp, report_specs, verbose=verbose)
        else:
            return self.get_relevant_enc_files(exp_id, exp, report_specs, reps, verbose=verbose)
            

    def file_size_string(self,size):
        '''Returns a file size string in G, M denomiations or else the number itself.'''
        if size >= 1000000000:
            return "%.2fG" % (size / 1000000000.0)
        if size >= 1000000:
            return "%.2fM" % (size / 1000000.0)
        if size >= 100000:
            return "%.1fK" % (size / 1000.0)
        return  size


    def gather_job_tree(self,job_id,job_cache,exclude_files=[],found_job_ids=[],level=0,verbose=False):
        '''Returns all parent jobs for the current job_id.  Recursive by default.
           job_id: the job at which to start the search.
           job_cache: can contain more jobs than just this search
           exclude_files: list of file endings to exclude (e.g. [ '_index.tgz' ]). 
           found_job_ids: cumulative list of all job_ids found in previous levels of a recursive search
           level: -1 means no recursion, 0 means start recursion at this level.  
           Returns (updated job_cache, list job_ids found for just one (or all recursive) levels)
        '''
        verbose=False # Verbose can be overkill!
        parent_job_ids = []
        
        if job_id in job_cache:
            job = job_cache[job_id]
        else:
            try:
                job =  dxpy.api.job_describe(job_id)
                job_cache[job_id] = job
            except:
                print >> sys.stderr, "WARNING: Could not find job for %s." % job_id
                return (job_cache,parent_job_ids)
                        
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
            print >> sys.stderr, "%s- Found %d inputs." % (' ' * level,len(file_inputs))

        # For each file input, verify it is for a file and then look for an accession.
        for inp in file_inputs:
            if not type(inp) == dict:
                if verbose:
                    print >> sys.stderr, "WARNING: Unexpected input_type: " + str(type(inp))
                continue # not a file input
            dxlink = inp.get("$dnanexus_link")
            if dxlink == None:
                continue
            if not type(dxlink) == dict:
                inp_fid = dxlink
            else:
                inp_fid = dxlink.get("id")
            try:
                file_dx_obj = dx.description_from_fid(inp_fid)
                job_id = file_dx_obj["createdBy"]["job"]
            except:
                # This is not an error... fastqs, fastas and chrom.sizes should have no job_ids!
                #if verbose:
                #    print >> sys.stderr, "WARNING: can't find file or job from "+ inp_fid
                continue

            # Determine jobs that shouldn't be followed: indexing jobs
            exclude = False
            for ending in exclude_files:
                if file_dx_obj["name"].endswith(ending):
                    exclude = True
                    break
            if exclude:
                continue

            # If recursing then 2 steps can have inputs from the same job, so don't count it twice 
            if job_id in found_job_ids:  
                continue
            # One step can have 2 inputs from the same job, so don't count it twice 
            if job_id in parent_job_ids:
                continue
            parent_job_ids.append(job_id)
            if verbose: 
                print >> sys.stderr, "%s* Found: %s" % (' ' * level,job_id)
        
        # Now recurse:
        if level >= 0:
            this_level_ids = parent_job_ids
            found_job_ids.extend( parent_job_ids )
            for job_id in this_level_ids:
                (job_cache, new_ids) = self.gather_job_tree(job_id,job_cache,exclude_files,found_job_ids,(level+1),verbose=verbose)
                for new_id in new_ids:
                    if new_id not in parent_job_ids:
                        parent_job_ids.append( new_id )
                    if new_id not in found_job_ids:
                        found_job_ids.append( new_id )
                        
        if verbose: 
            print >> sys.stderr, "%s- Level: %d, found job tree: %s" % (' ' * level,level,str(found_job_ids) )
            
        return (job_cache,parent_job_ids)
            
            
    def get_dx_special_metric(self,metric_id,file_dx_obj,metric_key,metric_def,verbose=False):
        '''Returns a mocked up 'metric' object from specialized definitions in DX.'''
        #verbose=False
        
        if verbose:
            print >> sys.stderr, "Looking for special metric: %s" % ( metric_id )
            
        metric = {}
        if "jobs" not in self.obj_cache["exp"]:
            self.obj_cache["exp"]["jobs"] = {}
        try:
            job_id = file_dx_obj["createdBy"]["job"]
            job =  dxpy.api.job_describe(job_id)
            self.obj_cache["exp"]["jobs"][job_id] = job
        except:
            print >> sys.stderr, "ERROR: Could not find job for %s." % file_dx_obj['name']
            return None
        #print >> sys.stderr, "file_dx_obj for '"+metric_id+"':"
        #print >> sys.stderr, json.dumps(file_dx_obj,indent=4)
        #print >> sys.stderr, "job:"
        #print >> sys.stderr, json.dumps(job,indent=4)
        #sys.exit(0)
          
        # In case job tree needs to be gathered
        excude_files = [ '.tgz' ] # All Reference files
        # All the fastas, fastqs and chrom.sizes should be excluded by having no job that created them
        # excude = [ '.fastq.gz', '.fq.gz', '.fastq', '.fq', '.chrom.sizes','.fasta.gz','.fa.gz','.fasta','.fa' ]
        all_jobs_ids = []
         
        for col in metric_def["columns"]:
            if verbose:
                print >> sys.stderr, "- Special metric column: '%s'" % ( col )
            if col == "Cost":
                cost = job.get("totalPrice")
                if cost != None:
                    metric[col] = "$" + str(round(cost,2))
            elif col == "Total Cost":
                # Plan: (1) get job input files all the way back to first job (excluding indexes), (2) sum all jobs
                total_cost = 0
                if len(all_jobs_ids) == 0:
                    (self.obj_cache["exp"]["jobs"], all_jobs_ids) = self.gather_job_tree(job_id,self.obj_cache["exp"]["jobs"], \
                                                                                              excude_files,[],verbose=verbose) 
                    if verbose:
                        print >> sys.stderr, "- '%s' from %d jobs" % ( col, len(all_jobs_ids) )
                for a_job_id in all_jobs_ids:
                    a_job = self.obj_cache["exp"]["jobs"][a_job_id]
                    total_cost += a_job.get("totalPrice",0)
                if total_cost != 0:
                    metric[col] = "$" + str(round(total_cost,2))
            elif col == "Time":
                beg = job.get("startedRunning") 
                end = job.get("stoppedRunning") 
                if beg != None and end != None:
                    duration = (end/1000.0 - beg/1000.0)
                    metric[col] = duration
            elif col == "Job Time":
                # Plan: (1) get job input files all the way back to first job (excluding indexes), (2) sum all jobs
                total_duration = 0
                if len(all_jobs_ids) == 0:
                    (self.obj_cache["exp"]["jobs"], all_jobs_ids) = self.gather_job_tree(job_id,self.obj_cache["exp"]["jobs"], \
                                                                                                excude_files,[],verbose=verbose) 
                    if verbose:
                        print >> sys.stderr, "- '%s' from %d jobs" % ( col, len(all_jobs_ids) )
                for a_job_id in all_jobs_ids:
                    a_job = self.obj_cache["exp"]["jobs"][a_job_id]
                    beg = a_job.get("startedRunning") 
                    end = a_job.get("stoppedRunning") 
                    if beg != None and end != None:
                        total_duration += (end/1000.0 - beg/1000.0)
                if total_duration != 0:
                    metric[col] = total_duration
            elif col == "CPU Time":
                # Plan: (1) get job input files all the way back to first job (excluding indexes), (2) sum all jobs
                total_duration = 0
                if len(all_jobs_ids) == 0:
                    (self.obj_cache["exp"]["jobs"], all_jobs_ids) = self.gather_job_tree(job_id,self.obj_cache["exp"]["jobs"], \
                                                                                                excude_files,[],verbose=verbose) 
                    if verbose:
                        print >> sys.stderr, "- '%s' from %d jobs" % ( col, len(all_jobs_ids) )
                for a_job_id in all_jobs_ids:
                    a_job = self.obj_cache["exp"]["jobs"][a_job_id]
                    beg = a_job.get("startedRunning") 
                    end = a_job.get("stoppedRunning") 
                    if beg != None and end != None:
                        duration = (end/1000.0 - beg/1000.0)
                    try:
                        instance = a_job["instanceType"]
                        cpus = int(a_job["instanceType"].split('x')[1])
                    except:
                        cpus = 1
                    total_duration = total_duration + (duration * cpus) 
                if total_duration != 0:
                    metric[col] = total_duration
                
            elif col == "File Size":
                size = file_dx_obj.get("size")
                if size != None:
                    metric[col] = size
        if verbose:
            print >> sys.stderr, "Special metric '"+metric_id+"':"
            print >> sys.stderr, json.dumps(metric,indent=4)
        return metric
            

    def get_dx_metrics(self,exp_id,target_files,report_specs,verbose=False):
        '''Returns a list of tuples of metrics objects from DX for reporting statistics.'''
        #verbose=True

        metrics = []
        combo_metrics = []
        metric_ids = []
        self.obj_cache["exp"]["metrics"] = {}
        for (out_type,sub_type,rep_tech,fid) in target_files:
            file_dx_obj = self.obj_cache["exp"]["files"][fid]
            
            if sub_type != None:
                metric_keys = report_specs[sub_type]["metrics"]
            else:
                metric_keys = report_specs[out_type]["metrics"]
            
            assert metric_keys != None
            if not isinstance(metric_keys,list):
                metric_keys = [ metric_keys ]
            
            for metric_key in metric_keys:                
                metric_defs = self.METRIC_DEFS[metric_key]
                
                # Make only one 'combined' metric from 2 reps.
                if metric_defs["per"] == "experiment":
                    metric_id = exp_id + '/' + metric_key
                else: 
                    metric_id = fid + '/' + metric_key 
                if metric_id in metric_ids:  
                    continue
                 
                if "special_metric" in metric_defs and metric_defs["special_metric"]:
                    if verbose:
                        print >> sys.stderr, "- Getting special metric for %s %s %s" % (exp_id,rep_tech,metric_key)
                    metric = self.get_dx_special_metric(metric_id,file_dx_obj,metric_key,metric_defs,verbose=verbose)
                else:
                    if verbose:
                        print >> sys.stderr, "- Getting qc metric for %s %s %s" % (exp_id,rep_tech,metric_key)
                    dx_key = metric_defs.get("dx_key",metric_key)
                    metric = file_dx_obj['QC'].get(dx_key)
                    # TODO: Add conversion of DX to ENC qc_metric objects if necessary
                if metric == None:
                    if verbose:
                        print >> sys.stderr, "  - Experiment "+exp_id+" "+file_dx_obj['name']+" has no " + metric_key + \
                                                                                                " quality_metric objects."
                    continue
                metric_ids.append(metric_id)
                if metric_defs["per"] == "experiment":
                    if verbose:
                        print >> sys.stderr, "  * Experiment "+exp_id+" found 'combined' " + metric_key + " quality_metric."
                    combo_metrics.append((metric_key,"combined",metric_id))
                else:
                    if verbose:
                        print >> sys.stderr, "  * Experiment "+exp_id+" found "+rep_tech+" " + metric_key + " quality_metric."
                    metrics.append((metric_key,rep_tech,metric_id))
                self.obj_cache["exp"]["metrics"][metric_id] = metric

        if len(combo_metrics) > 0:
            metrics.extend(combo_metrics)

        if verbose:
            print >> sys.stderr, "Encoded metrics:"
            for (metric_key,rep_tech,metric_id) in metrics:
                print >> sys.stderr, "%s %s %s" % (metric_key,rep_tech,metric_id)
        return metrics

    def total_step_runs(self, file):
        # this will recursively return a list of step runs
        # rewrite this so that it will gather all step runs for an experiment or for just 1 replicate of an experiment
        # the set of step runs will be ones that are in files or in the qc metric in the files
        # notes I'm given a file and I can go from the file into the parent experiment
        derived_from = file.get("derived_from", [])
        step_runs = []
        for fi in derived_from:
            if type(fi) is unicode:
              temp = fi
            else:
              temp = fi["@id"]
            # get the derived from file
            file_obj = encd.lookup_json(temp)
            if file_obj.get("output_type", "") == "genome index":
                continue
            step_run = file_obj.get("step_run")
            if step_run is None or step_run in step_runs:
                # we don't want this one
                continue
            # this means that we found a step run that we have not yet found
            step_runs.append(step_run)
            step_runs += self.total_step_runs(fi)
        return step_runs

    def get_enc_special_metric(self,metric_id,file_enc_obj,metric_key,metric_def,verbose=False):
        '''Returns a mocked up 'metric' object from specialized definitions in encodeD.'''
        #verbose=True
        
        metric = {}
        notes = None
        file_notes = file_enc_obj.get("notes")
        if file_notes != None:
            try:
                notes = json.loads(file_notes)
            except:
                notes = None
            
        for col in metric_def["columns"]:
            if col == "Cost":
                if notes != None:
                    cost = notes.get("dx_cost")
                    if cost == None:
                        step_notes = file_enc_obj["step_run"].get("notes")
                        if step_notes != None:
                            try:
                                cost = json.loads(step_notes).get("dx_cost")
                            except:
                                cost = None
                    if cost != None:
                        metric[col] = cost
            elif col == "Total Cost":
                # first need all the step runs
                # gather all the step runs
                total = 0
                step_runs = self.total_step_runs(file_enc_obj)
                step_runs = list(set(step_runs))
                for run in step_runs:
                    step = encd.lookup_json(run)
                    step_notes = step.get("notes")
                    if step_notes is not None:
                        try:
                            cost = json.loads(step_notes).get("dx_cost")
                            cost = float(cost.lstrip("$"))
                        except:
                            cost = None
                        if cost is not None:
                            total += cost
                metric[col] = total
            elif col == "Time":
                step_run = encd.lookup_json(file_enc_obj["step_run"])
                if step_run != None:
                    step_details = step_run["dx_applet_details"][0]
                    beg_time = step_details.get("started_running")
                    end_time = step_details.get("stopped_running")
                    if beg_time != None and end_time != None:
                        beg_dt = datetime.strptime(beg_time, '%Y-%m-%dT%H:%M:%SZ')
                        end_dt = datetime.strptime(end_time, '%Y-%m-%dT%H:%M:%SZ')
                        duration = (end_dt - beg_dt)
                        if duration != None:
                            metric[col] = duration.total_seconds()
            elif col == "Job Time":
                step_runs = self.total_step_runs(file_enc_obj)
                step_runs = list(set(step_runs))
                total = 0
                for run in step_runs:
                    step_run = encd.lookup_json(file_enc_obj["step_run"])
                    if step_run != None:
                        step_details = step_run["dx_applet_details"][0]
                        beg_time = step_details.get("started_running")
                        end_time = step_details.get("stopped_running")
                        if beg_time != None and end_time != None:
                            beg_dt = datetime.strptime(beg_time, '%Y-%m-%dT%H:%M:%SZ')
                            end_dt = datetime.strptime(end_time, '%Y-%m-%dT%H:%M:%SZ')
                            duration = (end_dt - beg_dt)
                            if duration != None:
                                total += duration.total_seconds()
                metric[col] = total
            elif col == "CPU Time": # TODO:
                pass
                #print >> sys.stderr, "WARNING: '%s' not yet supported for '%s'" % (col,metric_id)
            elif col == "File Size":
                metric[col] = file_enc_obj["file_size"]
        if verbose:
            print >> sys.stderr, "Special metric '"+metric_id+"':"
            print >> sys.stderr, json.dumps(metric,indent=4)
        return metric
            

    def get_enc_metrics(self,exp_id,target_files,report_specs,verbose=False):
        '''Returns a list of tuples of metrics objects from encodeD for reporting statistics.'''
        #verbose=True
        
        metrics = []
        combo_metrics = []
        metric_ids = []
        self.obj_cache["exp"]["metrics"] = {}
        for (out_type,sub_type,rep_tech,acc) in target_files:
            file_enc_obj = self.obj_cache["exp"]["files"][acc]
            if sub_type != None:
                metric_keys = report_specs[sub_type]["metrics"]
            else:
                metric_keys = report_specs[out_type]["metrics"]
            
            assert metric_keys != None
            if not isinstance(metric_keys,list):
                metric_keys = [ metric_keys ]
            
            for metric_key in metric_keys:                
                metric_defs = self.METRIC_DEFS[metric_key]
                if "special_metric" in metric_defs and metric_defs["special_metric"]:
                    if verbose:
                        print >> sys.stderr, "- Getting special metric for %s %s %s" % (exp_id,rep_tech,metric_key)
                    metric_id = acc + '/' + metric_key 
                    if metric_id not in metric_ids:  # Note that this will result in only one combined metric from 2 reps.
                        metric_ids.append(metric_id) 
                        metric = self.get_enc_special_metric(metric_id,file_enc_obj,metric_key,metric_defs,verbose=verbose)
                        if metric != None:
                            if metric_defs["per"] == "experiment":
                                combo_metrics.append((metric_key,"combined",metric_id))
                            else:
                                metrics.append((metric_key,rep_tech,metric_id))
                            self.obj_cache["exp"]["metrics"][metric_id] = metric
                else:
                    if verbose:
                        print >> sys.stderr, "- Getting qc metric for %s %s %s" % (exp_id,rep_tech,metric_key)
                    file_qc_metrics = file_enc_obj.get('quality_metrics')
                    if file_qc_metrics == None or len(file_qc_metrics) == 0:
                        if verbose:
                            print >> sys.stderr, "  - Experiment "+exp_id+" file "+acc+" has no quality_metric objects."
                        continue
                    #elif verbose:
                    #    print >> sys.stderr, json.dumps(file_qc_metrics,indent=4,sort_keys=True)
                    qc_type = string.capwords(metric_defs.get("name",metric_key),'_').replace('_','') + "QualityMetric"
                    for metric in file_qc_metrics:
                        #print >> sys.stderr, json.dumps(metric,indent=4)
                        if verbose:
                            print >> sys.stderr, "  examining %s" % (metric["@type"])
                        if qc_type not in metric["@type"]:
                            if verbose:
                                print >> sys.stderr, "  - %s not found in %s." % (qc_type,metric["@type"])
                            continue
                        metric_id = metric["@id"]
                        if metric_id not in metric_ids:  # Note that this will result in only one combined metric from 2 reps.
                            if verbose:
                                print >> sys.stderr, "  adding id %s" % (metric_id)
                            metric_ids.append(metric_id) 
                            if metric_defs["per"] == "experiment":
                                if verbose:
                                    print >> sys.stderr, "  combining %s" % (metric["@type"])
                                combo_metrics.append((metric_key,"combined",metric_id))
                            else:
                                if verbose:
                                    print >> sys.stderr, "  not combining %s" % (metric["@type"])
                                metrics.append((metric_key,rep_tech,metric_id))
                            self.obj_cache["exp"]["metrics"][metric_id] = metric
        if len(combo_metrics) > 0:
            metrics.extend(combo_metrics)

        if verbose:
            print >> sys.stderr, "Encoded metrics:"
            for (metric_key,rep_tech,metric_id) in metrics:
                print >> sys.stderr, "%s %s %s" % (metric_key,rep_tech,metric_id)
        return metrics


    def get_metrics(self,exp_id,target_files,report_specs,verbose=False):
        '''Returns a list of tuples of metrics objects for reporting statistics.'''
        if self.data_mine == "DX":
            return self.get_dx_metrics(exp_id,target_files,report_specs,verbose=verbose)
        else:
            return self.get_enc_metrics(exp_id,target_files,report_specs,verbose=verbose)
        

    def type_headers(self,section_type,report_specs):
        '''Generates the headers for a type section from report specs.'''
        type_header = ""
        metric_keys = report_specs[section_type]["metrics"]
        if not isinstance(metric_keys,list):
            metric_keys = [ metric_keys ]
        for metric_key in metric_keys:
            metric_defs = self.METRIC_DEFS[metric_key]
            for col in metric_defs['columns']:
                if "headings" in metric_defs and col in metric_defs["headings"].keys():
                    heading = metric_defs["headings"][col]
                else:
                    heading = col
                if "prepend" in metric_defs:
                    if metric_defs["prepend"] == "type":
                        heading = section_type.capitalize() + ' ' + heading
                # Note: could try fancy stuff to stack column headers ( header vs "head\ner").
                #if len(heading) > 10 and len(heading.split()) > 1:
                #    new_head = '"'
                #    cur_len = 0
                #    for word in heading.split():
                #        if (cur_len + len(word)) > 10 and cur_len > 0:
                #            new_head += '\n' + word
                #            cur_len = len(word)
                #        elif cur_len > 0:
                #            new_head += ' ' + word
                #            cur_len += 1 + len(word)
                #        else:
                #            new_head += word
                #            cur_len = len(word)
                #    heading = new_head + '"'
                    
                type_header += '\t' + heading
            
        return type_header

    def print_headers(self,report_specs):
        '''Prints headers from report specs.'''
        Header = "# Experiment\treplicate"
        for out_type in report_specs["output_types"]:
            types_list = [ out_type ]
            if "subtypes" in report_specs[out_type]:
                types_list = report_specs[out_type]["subtypes"]
            for sub_type in types_list:
                Header += self.type_headers(sub_type,report_specs)
            
        print Header

    def init_type_totals(self,section_type,report_specs):
        '''Intitializes totals for a type section from report specs.'''
        metric_keys = report_specs[section_type]["metrics"]
        if not isinstance(metric_keys,list):
            metric_keys = [ metric_keys ]
        for metric_key in metric_keys:
            metric_defs = self.METRIC_DEFS[metric_key]
            if self.exp_starts == -1 and metric_defs["per"] == "experiment":
                self.exp_starts = len(self.min_stats)
            for col in metric_defs['columns']:
                self.min_stats.append(999999999) # TODO: deal with non-numeric columns
                self.max_stats.append(-999999999)
                self.cum_stats.append(0.0)
                self.col_counts.append(0)
                if col in metric_defs.get('sumable',[]):
                    self.sum_these_cols.append(True)
                else:
                    self.sum_these_cols.append(False)
        
    def init_totals(self,report_specs):
        '''Intitializes totals from report specs.'''
        self.min_stats = []
        self.max_stats = []
        self.cum_stats = []
        self.col_counts = []
        self.col_formats = {}
        self.rep_count = 0
        self.exp_count = 0
        self.exp_starts = -1
        self.not_tabulated_cols = [] # Some columns are not numeric
        self.sum_these_cols = []     # Some columns can be summed
        for out_type in report_specs["output_types"]:
            types_list = [ out_type ]
            if "subtypes" in report_specs[out_type]:
                types_list = report_specs[out_type]["subtypes"]
            for sub_type in types_list:
                self.init_type_totals(sub_type,report_specs)

    def format_value(self,val,format_type):
        '''Returns formatted string.'''
        if format_type == "file_size":
            return self.file_size_string(val)
        elif format_type == "duration":
            return dx.duration_string(val,include_seconds=False)
        elif format_type == "percent":
            if isinstance(val,int):
                return str(val)+"%"
            else:
                return "%.1f%%" % val
        elif format_type == "dollar":
            if isinstance(val,int):
                return "$"+str(val)
            else:
                return "$"+ "%.2f" % val
        else:
            print >> sys.stderr, "Unknown format type: '%s'" % (format_type)
            return str(val) 
            
            
    def print_metrics(self,exp_id,target_files,metrics,report_specs, verbose=False):
        '''Prints the requested metrics.'''
        #verbose=True
        # exp rep1_1:   n n n _ _ _
        #     rep2_1:   n n n _ _ _
        #     combined: _ _ _ n n n
        # totals min:   n n n n n n
        #        max:   n n n n n n
        #        mean:  n n n n n n
        Line = exp_id
        last_rep_tech = None
        for (metric_key,rep_tech,metric_id) in metrics:
            metric_defs = self.METRIC_DEFS[metric_key]
            
            # starting new line?
            new_line = False
            if last_rep_tech == None:
                new_line = True
            elif rep_tech != last_rep_tech:
                print Line
                Line = '' # Successive lines do not repeat exp_id
                new_line = True
            if new_line:
                last_rep_tech = rep_tech
                start_col = 0
                
            # What kind of new line?
            if new_line:
                if metric_defs["per"] == "replicate":
                    self.rep_count += 1
                    Line += '\t'+rep_tech+':'
                if metric_defs["per"] == "experiment":
                    start_col = self.exp_starts
                    Line += '\tcombined:'
                for col in range(start_col):
                    Line += '\t'
                    
            metric = self.obj_cache["exp"]["metrics"][metric_id]
            cur_col = start_col
            for col in metric_defs['columns']:
                val = metric.get(col)
                if val == None:
                    if verbose:
                        print >> sys.stderr, "EXP: %s col: %d is empty" % (exp_id,cur_col) 
                    Line += '\t'
                elif "format" in metric_defs and col in metric_defs["format"]:
                    format_type = metric_defs["format"][col]
                    self.col_formats[cur_col] = format_type
                    Line += '\t'+self.format_value(val,format_type)
                elif isinstance(val,str) or isinstance(val,unicode):
                    Line += '\t'+val
                    if verbose:
                        print >> sys.stderr, "EXP: %s col: %d is str val: '%s'" % (exp_id,cur_col,val) 
                    # can convert to number?
                    if val.startswith('%') or val.endswith('%'):
                        val = val.strip('%')
                        self.col_formats[cur_col] = "percent"
                    elif val.startswith('$'):
                        val = val.strip('$')
                        self.col_formats[cur_col] = "dollar"
                    else:
                        val = None
                    if val != None:
                        try:
                            val = int(val)
                        except:
                            try:
                                val = float(val)
                            except:
                                if verbose:
                                    print >> sys.stderr, "Couldn't convert '"+str(val)+"' to number."
                                val = None # No totaling
                    if val == None:
                        if cur_col not in self.not_tabulated_cols:
                            self.not_tabulated_cols.append(cur_col)
                elif isinstance(val,int) or isinstance(val,float):
                    Line += '\t'+str(val)
                else: # TODO: ignoring for now
                    if verbose:
                        print >> sys.stderr, "EXP: %s col: %d is not str, int or float but '%s'" % (exp_id,cur_col,type(val))
                    if cur_col not in self.not_tabulated_cols:
                        self.not_tabulated_cols.append(cur_col)
                    val = None
                if val != None:
                    if  self.min_stats[cur_col] > val or self.min_stats[cur_col] == 999999999:
                        self.min_stats[cur_col] = val
                    if  self.max_stats[cur_col] < val or self.max_stats[cur_col] == -999999999:
                        self.max_stats[cur_col] = val
                    self.cum_stats[cur_col] += val
                    self.col_counts[cur_col] += 1
                        
                cur_col += 1  # Done with column
                start_col = cur_col
        if Line != "":
            print Line
        self.exp_count += 1 # Done with experiment
                    
                                 
    def print_totals(self):
        '''Prints accumulated totals.'''

        # min:
        Line = "# Totals\tmin:"
        cur_col = 0
        for val in self.min_stats:
            if cur_col in self.not_tabulated_cols or val == 999999999 or self.col_counts[cur_col] == 0:
                Line += '\t'
            else:
                if cur_col in self.col_formats.keys():
                    Line += '\t'+self.format_value(val,self.col_formats[cur_col])
                else:                 
                    Line += '\t'+str(val)
            cur_col += 1
        print Line

        # max:
        Line = "# \tmax:"
        cur_col = 0
        for val in self.max_stats:
            if cur_col in self.not_tabulated_cols or val == -999999999 or self.col_counts[cur_col] == 0:
                Line += '\t'
            else:
                if cur_col in self.col_formats.keys():
                    Line += '\t'+self.format_value(val,self.col_formats[cur_col])
                else:                 
                    Line += '\t'+str(val)
            cur_col += 1
        print Line
        
        # mean:
        Line = "# \tmean:"
        cur_col = 0
        #print >> sys.stderr, "exp_starts: %d" % self.exp_starts
        for val in self.cum_stats:
            if cur_col in self.not_tabulated_cols or self.col_counts[cur_col] == 0:
                Line += '\t'
            else:
                if self.col_counts[cur_col] > 1:
                    val = val / self.col_counts[cur_col]
                #if self.exp_starts == -1 or cur_col < self.exp_starts:
                #    if self.rep_count > 0:
                #        val = val / self.rep_count # TODO: make sure it is a float!
                #else: 
                #    if self.exp_count > 0:
                #        val = val / self.exp_count
                        
                if cur_col in self.col_formats.keys():
                    Line += '\t'+self.format_value(val,self.col_formats[cur_col])
                elif val > 9999:
                    Line += '\t%.1f' % val
                elif val > 999:
                    Line += '\t%.2f' % val
                elif val > 99:
                    Line += '\t%.3f' % val
                elif val > 9:
                    Line += '\t%.4f' % val
                else:
                    Line += '\t%f' % val
            cur_col += 1
        print Line

        
        # Sum:
        Line = "# \tsum:"
        cur_col = 0
        #print >> sys.stderr, "exp_starts: %d" % self.exp_starts
        print_sums = False
        for val in self.cum_stats:
            if not self.sum_these_cols[cur_col] or cur_col in self.not_tabulated_cols or self.col_counts[cur_col] == 0:
                Line += '\t'
            else:
                if self.col_counts[cur_col] > 1:
                    val = val
                    print_sums = True
                #if self.exp_starts == -1 or cur_col < self.exp_starts:
                #    if self.rep_count > 0:
                #        val = val / self.rep_count # TODO: make sure it is a float!
                #else: 
                #    if self.exp_count > 0:
                #        val = val / self.exp_count
                        
                if cur_col in self.col_formats.keys():
                    Line += '\t'+self.format_value(val,self.col_formats[cur_col])
                elif val > 9999:
                    Line += '\t%.1f' % val
                elif val > 999:
                    Line += '\t%.2f' % val
                elif val > 99:
                    Line += '\t%.3f' % val
                elif val > 9:
                    Line += '\t%.4f' % val
                else:
                    Line += '\t%f' % val
            cur_col += 1
        if print_sums:
            print Line
            
        # counts:
        print "# Experiments:\t" + str(self.exp_count)
        print "# Replicates:\t" + str(self.rep_count)


    def run(self):
        '''Runs mission_log from start to finish using command line arguments.'''
        args = self.get_args()
            
        self.server_key = args.server
        encd.set_server_key(self.server_key) # TODO: change to self.encd = Encd(self.server_key)
        #(self.authid, self.authpw, self.server) = encd.find_keys(self.server_key)
        self.data_mine = "encodeD"
        if args.dx:
            self.data_mine = "DX"
            self.proj_name = dx.env_get_current_project()
            self.project = dx.get_project(self.proj_name)
            self.folder = args.folder
            print >> sys.stderr, "Using %s as data mine" % self.data_mine
        
        # Look up report type
        if args.report_type not in self.REPORT_SPECS.keys():
            print >> sys.stderr, "Report type %s is not supported" % (args.report_type)
            sys.exit(1)
        if self.data_mine not in self.REPORT_SPECS[args.report_type]["sources"]:
            print >> sys.stderr, "ERROR: Source '%s' is not supported for '%s'." % (self.data_mine, args.report_type)
            sys.exit(1)

        if args.genome != None:
            self.genome = args.genome
        else:
            self.genome_warning = True
            
        #if args.verbose:
        print >> sys.stderr, "== Running mission_log from [%s] server ==" % (self.server_key)
        
        self.exp_ids = self.load_exp_list(args.experiments,args.file,verbose=args.verbose)
        if len(self.exp_ids) == 0:
            print >> sys.stderr, "No experiment id's requested."
            self.ap.print_help()
            sys.exit(1)       

        exp_count = 0
        for exp_id in self.exp_ids:
            sys.stdout.flush() # Slow running job should flush to piped log
            self.exp_id = exp_id
            self.obj_cache["exp"] = {}  # clear exp cache, which will hold exp specific wf_run and step_run objects
            # Lookup experiment type from encoded, based on accession
            print >> sys.stderr, "Working on %s..." % self.exp_id
            self.exp = encd.get_exp(self.exp_id,must_find=True)
            if self.exp == None or self.exp["status"] == "error":
                print >> sys.stderr, "Unable to locate experiment %s in encoded (%s)" % (self.exp_id, self.server_key)
                continue
                
            # Type of experiment?
            exp_type = encd.get_assay_type(self.exp_id,self.exp)
            if self.exp_type != None and self.exp_type != exp_type:
                print >> sys.stderr, "WARNING: Experiment type '%s' does not match previous type(s) '%s'.  Skipping!" % \
                                                                                            (exp_type,self.exp_type)
                continue
            if self.exp_type == None:
                self.exp_type = exp_type
                if self.exp_type not in self.EXPERIMENT_TYPES_SUPPORTED:
                    print >> sys.stderr, "ERROR: Experiment %s has unsupported assay type of '%s'" % (self.exp_id,self.exp_type)
                    sys.exit(1)

                if self.exp_type not in self.REPORT_SPECS[args.report_type]:
                    print >> sys.stderr, "ERROR: Report type %s is not supported for experiment type '%s'" % \
                                                                                    (args.report_type,self.exp_type)
                    sys.exit(1)
                self.report_specs = self.REPORT_SPECS[args.report_type][self.exp_type]
                    
                # We have enough to begin our report so...
                self.print_headers(self.report_specs)
                self.init_totals(self.report_specs)
                
            # Now for each experiment:
            
            # Get reps from encoded
            self.reps = self.find_replicates(self.exp_id, self.exp)
                
            # Get files from encoded
            target_files = self.get_relevant_files(self.exp_id, self.exp, self.report_specs, self.reps,verbose=args.verbose)
            
            if target_files == None or len(target_files) == 0:
                print >> sys.stderr, "WARNING: No files available on encodeD for %s" % self.exp_id
                continue
            
            # Get metrics per file
            metrics = self.get_metrics(self.exp_id,target_files,self.report_specs,verbose=args.verbose)
            if metrics == None or len(metrics) == 0:
                print >> sys.stderr, "WARNING: No metrics available on %s for %s" % (self.data_mine,self.exp_id)
                continue
            
            # Print metrics
            self.print_metrics(self.exp_id,target_files,metrics,self.report_specs)
            
        # print any totals available
        if self.report_specs:     
            self.print_totals()
        print >> sys.stderr, "(finished)"


if __name__ == '__main__':
    '''Run from the command line.'''
    mission_log = Mission_log()
    mission_log.run()

