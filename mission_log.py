#!/usr/bin/env python2.7
# mission_log.py 0.1.0
#
# Mission_log is meant to report on experiments that have been posted to encoded.
# It expects a list of experiment accessions to report on.
#
# 1) List of experiments
# 2) each experiment
#   a) find star bams and quants 
#      dxencode.get_enc_exp_files(exp_obj,["transcriptome alignments","gene quantifications"])
#      - order by replicate!
#   b) find qc metric associated:
#      splashdown.enc_lookup_json("/star_quality_metric/dghdf",frame='object',must_find=False)
#   c) read qc metric for each bam
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
import commands

import dxpy
import dxencode

class Mission_log(object):
    '''
    Mission_log module reports information for a set of experiments posted to ENCODEd.
    '''
    TOOL_IS = 'mission_log'
    HELP_BANNER = "Reports statistics for a set of experiments posted to ENCODEd. " 
    ''' This help banner is displayed by get_args.'''

    SERVER_DEFAULT = 'www'
    '''This the default server to report from.'''

    EXPERIMENT_TYPES_SUPPORTED = [ 'long-rna-seq', 'small-rna-seq'] #, 'rampage','dnase','dna-me','chip-seq' ]
    '''This module supports only these experiment (pipeline) types.'''

    REPORT_DEFAULT = 'star-mad'
    '''This the default report type.'''
    REPORTS_SUPPORTED = [ 'star-mad' ]
    '''This module supports only these report types.'''

    ASSEMBLIES_SUPPORTED = { "hg19": "hg19", "hg38": "GRCh37", "mm10": "mm10" }
    '''This module supports only these assemblies.'''

    ANNOTATIONS_SUPPORTED = [ 'V19', 'M2', 'M3', 'M4' ]
    '''This module supports only these annotations.'''

    REPORT_SPECS = {
         'star-mad': {
             "long-rna-seq":  { "output_types": [ "transcriptome alignments",        "gene quantifications" ],
                                "qc_metrics":   { "transcriptome alignments": 'star',"gene quantifications": 'mad' } },
             "small-rna-seq": { "output_types": [ "alignments",        "gene quantifications" ],
                                "qc_metrics":   { "alignments": 'star',"gene quantifications": 'mad' } }
         }
    }
    '''For each report type, these are the mappings for file 'output_types' to 'quality_metric' types.'''
    
    METRIC_DEFS = {
       'star': {"per": "replicate",
                "columns": [
                            "Number of input reads",
                            "Uniquely mapped reads number",
                            "Uniquely mapped reads %" ],
                "headings": {
                                "Number of input reads":        "Reads",
                                "Uniquely mapped reads number": "Unique Reads",
                                "Uniquely mapped reads %":      "Unique Read %"}
              },
       'mad': { "per": "experiment",
                "columns": [
                            "MAD of log ratios",
                            "Pearson correlation",
                            "Spearman correlation" ],
              }
    }

    def __init__(self):
        '''
        Splashdown expects one or more experiment ids as arguments and will find, document
        and post files in the associated directory.
        '''
        self.args = {} # run time arguments
        self.server_key = 'www'
        self.exp = {}  # Will hold the encoded exp json
        self.exp_id = None
        self.exp_type = None  # Will hold the experiment's assay_type, normalized to known tokens.
        self.genome = None  # TODO: need way to determine genome before any posts occur!
        self.annotation = None  # TODO: if appropriate, need way to determine annotation
        self.report_specs = {} # Unless the proper report specifications are discovered, a report can't be run
        self.obj_cache = {} # certain things take time to find or create and are needed multiple times
        print >> sys.stderr, " "


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
                        required=False)

        ap.add_argument('-f','--file',
                        help="A file containing experiment accession (first word in line, # comments ignored)",
                        default=None,
                        required=False)

        ap.add_argument('-r','--report-type',
                        help="The report type to print (default: '" + self.REPORT_DEFAULT + "')",
                        default=self.REPORT_DEFAULT,
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
            self.obj_cache["exp"]["full_mapping"] = dxencode.get_full_mapping(exp_id,exp,key=self.server_key)
        replicates = dxencode.get_reps_from_enc(exp_id, exp=exp, full_mapping=self.obj_cache["exp"]["full_mapping"], \
                                                                                                    key=self.server_key)
        for rep in replicates:
            if self.genome == None:
                if rep['organism'] in dxencode.GENOME_DEFAULTS:
                    self.genome = dxencode.GENOME_DEFAULTS[rep['organism']]
                else:
                    print >> sys.stderr, "Organism %s not currently supported" % rep['organism']
                    sys.exit(1)
            elif self.genome != dxencode.GENOME_DEFAULTS[rep['organism']]:
                print >> sys.stderr, "Mixing genomes in one mission_log run not supported %s and %s" % \
                                                    (self.genome, dxencode.GENOME_DEFAULTS[rep['organism']])
                sys.exit(1)
                
        if verbose:
            print >> sys.stderr, "Replicates:"
            print >> sys.stderr, json.dumps(replicates,indent=4)
        return replicates


    def get_relevant_files(self, exp_id, exp, report_specs, reps, verbose=False):
        '''Returns list of enc file_objects for selected files available on encoded.'''
        #verbose=True
        enc_files = []
        self.obj_cache["exp"]["files"] = {}
        
        files = dxencode.get_enc_exp_files(exp,report_specs["output_types"],key=self.server_key)
        if len(files) == 0:
            if verbose:
                print >> sys.stderr, "Experiment "+exp_id+" has no results to report."
            return files

        # How to place them in order?
        for out_type in report_specs["output_types"]:
            for rep in reps:
                for file_obj in files:
                    if file_obj.get('output_type') != out_type:
                        continue
                    br = file_obj['replicate']['biological_replicate_number']
                    tr = file_obj['replicate']['technical_replicate_number']
                    if br != rep['br'] or tr != rep['tr']:
                        continue
                    rep_tech = "rep%d_%d" % (br,tr)
                    acc = file_obj['accession'] 
                    enc_files.append((out_type,rep_tech,acc))
                    self.obj_cache["exp"]["files"][acc] = file_obj

        if verbose:
            print >> sys.stderr, "Encoded files:"
            for (out_type,rep_tech,acc) in enc_files:
                print >> sys.stderr, "'%s' %s %s %s" % \
                        (out_type,rep_tech,acc,self.obj_cache["exp"]["files"][acc]['submitted_file_name'])
        return enc_files


    def get_qc_metrics(self,exp_id,enc_files,report_specs,verbose=False):
        '''Returns a list of tuples of qc metrics for reporting statistics.'''
        #verbose=True
        qc_metrics = []
        qc_ids = []
        self.obj_cache["exp"]["qc_metrics"] = {}
        for (out_type,rep_tech,acc) in enc_files:
            file_obj = self.obj_cache["exp"]["files"][acc]
            file_qc_metrics = file_obj.get('quality_metrics')
            if file_qc_metrics == None or len(file_qc_metrics) == 0:
                if verbose:
                    print >> sys.stderr, "Experiment "+exp_id+" file "+acc+" has no quality_metric objects."
                continue
            qc_key = report_specs["qc_metrics"].get(out_type)
            assert qc_key != None
            for qc_metric in file_qc_metrics:
                #print >> sys.stderr, json.dumps(qc_metric,indent=4)
                if qc_key+"_quality_metric" not in qc_metric["@type"]:
                    continue
                qc_id = qc_metric["@id"]
                if qc_id not in qc_ids:  # Note that this will automatically give only one combined metric from 2 reps.
                    qc_ids.append(qc_id) 
                    qc_metrics.append((qc_key,rep_tech,qc_id))
                    self.obj_cache["exp"]["qc_metrics"][qc_id] = qc_metric

        if verbose:
            print >> sys.stderr, "Encoded qc_metrics:"
            for (qc_key,rep_tech,qc_id) in qc_metrics:
                print >> sys.stderr, "%s %s %s" % (qc_key,rep_tech,qc_id)
        return qc_metrics


    def print_header(self,report_specs):
        '''Generates a header from report specs.'''
        Header = "# Experiment\treplicate"
        # TODO: fancy stuf to stack column headers
        for out_type in report_specs["output_types"]:
            qc_key = report_specs["qc_metrics"][out_type]
            metric_defs = self.METRIC_DEFS[qc_key]
            for col in metric_defs['columns']:
                if "headings" in metric_defs and col in metric_defs["headings"].keys():
                    heading = metric_defs["headings"][col]
                else:
                    heading = col
                Header += '\t' + heading
        print Header

    def init_totals(self,report_specs):
        '''Intitializes totals from report specs.'''
        self.min_stats = []
        self.max_stats = []
        self.cum_stats = []
        self.rep_count = 0
        self.exp_count = 0
        self.exp_starts = -1
        self.not_tabulated_cols = [] # Some columns are not numeric
        for out_type in report_specs["output_types"]:
            qc_key = report_specs["qc_metrics"][out_type]
            metric_defs = self.METRIC_DEFS[qc_key]
            if self.exp_starts == -1 and metric_defs["per"] == "experiment":
                self.exp_starts = len(self.min_stats)
            for col in metric_defs['columns']:
                self.min_stats.append(999999999) # TODO: deal with non-numeric columns
                self.max_stats.append(-999999999)
                self.cum_stats.append(0.0)

    def print_metrics(self,exp_id,enc_files,qc_metrics,report_specs, verbose=False):
        '''Prints the requested metrics.'''
        #verbose=True
        # exp rep1_1:   n n n _ _ _
        #     rep2_1:   n n n _ _ _
        #     combined: _ _ _ n n n
        # totals min:   n n n n n n
        #        max:   n n n n n n
        #        mean:  n n n n n n
        Line = exp_id
        for (qc_key,rep_tech,qc_id) in qc_metrics:
            metric_defs = self.METRIC_DEFS[qc_key]
            start_col = 0
            if metric_defs["per"] == "replicate":
                Line += '\t'+rep_tech+':'
                self.rep_count += 1
            if metric_defs["per"] == "experiment":
                start_col = self.exp_starts
                Line += '\tcombined:'
                for col in range(start_col):
                    Line += '\t'
            qc_metric = self.obj_cache["exp"]["qc_metrics"][qc_id]
            cur_col = start_col
            for col in metric_defs['columns']:
                val = qc_metric.get(col)
                if val == None:
                    if verbose:
                        print >> sys.stderr, "EXP: %s col: %d is empty" % (exp_id,cur_col) 
                    Line += '\t'
                elif isinstance(val,str) or  isinstance(val,unicode):
                    Line += '\t'+val
                    if verbose:
                        print >> sys.stderr, "EXP: %s col: %d is str" % (exp_id,cur_col) 
                    # can convert to number?
                    if val.startswith('%') or val.endswith('%'):
                        val = val.strip('%')
                        try:
                            val = int(val)
                        except:
                            try:
                                val = float(val)
                            except:
                                if cur_col not in self.not_tabulated_cols:
                                    self.not_tabulated_cols.append(cur_col)
                                val = None
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
                cur_col += 1  # Done with column
            print Line
            Line = '' # Successive lines do not repeat exp_id
        self.exp_count += 1 # Done with experiment
                    
                                 
    def print_totals(self):
        '''Prints accumulated totals.'''

        # min:
        Line = "# Totals\tmin:"
        cur_col = 0
        for val in self.min_stats:
            if cur_col in self.not_tabulated_cols or val == 999999999:
                Line += '\t'
            else:
                Line += '\t'+str(val)
            cur_col += 1
        print Line

        # max:
        Line = "# \tmax:"
        cur_col = 0
        for val in self.max_stats:
            if cur_col in self.not_tabulated_cols or val == -999999999:
                Line += '\t'
            else:
                Line += '\t'+str(val)
            cur_col += 1
        print Line
        
        # mean:
        Line = "# \tmean:"
        cur_col = 0
        #print >> sys.stderr, "exp_starts: %d" % self.exp_starts
        for val in self.cum_stats:
            if cur_col in self.not_tabulated_cols:
                Line += '\t'
            else:
                if self.exp_starts != -1 and cur_col < self.exp_starts:
                    if self.rep_count > 0:
                        val = val / self.rep_count # TODO: make sure it is a float!
                else: 
                    if self.exp_count > 0:
                        val = val / self.exp_count
                if val > 9999:
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

        # counts:
        print "# Experiments:\t" + str(self.exp_count)
        print "# Replicates:\t" + str(self.rep_count)


    def run(self):
        '''Runs mission_log from start to finish using command line arguments.'''
        args = self.get_args()
            
        self.server_key = args.server
        self.authid, self.authpw, self.server = dxencode.processkey(self.server_key)
        
        if args.verbose:
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
            # 1) Lookup experiment type from encoded, based on accession
            if args.verbose:
                print >> sys.stderr, "Working on %s..." % self.exp_id
            self.exp = dxencode.get_exp(self.exp_id,must_find=True,key=self.server_key)
            if self.exp == None or self.exp["status"] == "error":
                print >> sys.stderr, "Unable to locate experiment %s in encoded (%s)" % (self.exp_id, self.server_key)
                continue
                
            # Type of experiment?
            exp_type = dxencode.get_assay_type(self.exp_id,self.exp)
            if self.exp_type != None and self.exp_type != exp_type:
                print >> sys.stderr, "Experiment type '%s' does not match previous type(s) '%s'.  Skipping!" % \
                                                                                            (exp_type,self.exp_type)
                continue
            if self.exp_type == None:
                self.exp_type = exp_type
                if self.exp_type not in self.EXPERIMENT_TYPES_SUPPORTED:
                    print >> sys.stderr, "Experiment %s has unsupported assay type of '%s'" % (self.exp_id,self.exp_type)
                    sys.exit(1)

                # Look up report type
                if args.report_type not in self.REPORTS_SUPPORTED or args.report_type not in self.REPORT_SPECS.keys():
                    print >> sys.stderr, "Report type %s is not supported" % (args.report_type)
                    sys.exit(1)
                if self.exp_type not in self.REPORT_SPECS[args.report_type]:
                    print >> sys.stderr, "Report type %s is not supported for experiment type '%s'" % \
                                                                                    (args.report_type,self.exp_type)
                    sys.exit(1)
                self.report_specs = self.REPORT_SPECS[args.report_type][self.exp_type]
                
                # We have enough to begin our report so...
                self.print_header(self.report_specs)
                self.init_totals(self.report_specs)
                
            # Now for each experiment:
            
            # Get files from encoded
            self.reps = self.find_replicates(self.exp_id, self.exp)
            enc_files = self.get_relevant_files(self.exp_id, self.exp, self.report_specs, self.reps,verbose=args.verbose)
            
            # Get qc_metrics per file
            qc_metrics = self.get_qc_metrics(self.exp_id,enc_files,self.report_specs,verbose=args.verbose)
            
            # Print metrics
            self.print_metrics(self.exp_id,enc_files,qc_metrics,self.report_specs)
            
        # print any totals available
        if self.report_specs:     
            self.print_totals()
        if args.verbose:
            print "(finished)"


if __name__ == '__main__':
    '''Run from the command line.'''
    mission_log = Mission_log()
    mission_log.run()

