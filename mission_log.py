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

    EXPERIMENT_TYPES_SUPPORTED = [ 'long-rna-seq', 'small-rna-seq', 'rampage'] #, 'rampage','dnase','dna-me','chip-seq' ]
    '''This module supports only these experiment (pipeline) types.'''

    REPORT_DEFAULT = 'star-mad'
    '''This the default report type.'''
    REPORTS_SUPPORTED = [ 'star-mad', 'tophat_times' ]
    '''This module supports only these report types.'''

    ASSEMBLIES_SUPPORTED = { "hg19": "hg19", "hg38": "GRCh37", "mm10": "mm10" }
    '''This module supports only these assemblies.'''

    ANNOTATIONS_SUPPORTED = [ 'V19', 'M2', 'M3', 'M4' ]
    '''This module supports only these annotations.'''

    REPORT_SPECS = {
         'star-mad': {
             "long-rna-seq":  { "output_types": [ "transcriptome alignments",        "gene quantifications" ],
                                "metrics":      { "transcriptome alignments": 'star',"gene quantifications": 'mad' } },
             "small-rna-seq": { "output_types": [ "alignments",        "gene quantifications" ],
                                "metrics":      { "alignments": 'star',"gene quantifications": 'mad' } },
             "rampage":       { "output_types": [ "alignments",        "gene quantifications" ],
                                "metrics":      { "alignments": 'star',"gene quantifications": 'mad' } }
         },
         'tophat_times': {
             "long-rna-seq":  { "output_types": [ "alignments" ],
                                "ending":       [ ("alignments", "star", "_star_genome.bam"), 
                                                  ("alignments", "tophat", "_tophat.bam"   ) ],
                                "metrics":      { "star": ["star2", "file_cost"], "tophat": "file_cost", } },
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
       'star2': {"per": "replicate",
                "name": "star",
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
    }

    def __init__(self):
        '''
        Mission_log generates reports from pipeline results that include quality_metrics that have already been 
        uploaded to encodeD.  
        '''
        self.args = {} # run time arguments
        self.server_key = 'www'
        self.exp = {}  # Will hold the encoded exp json
        self.exp_id = None
        self.exp_type = None  # Will hold the experiment's assay_type, normalized to known tokens.
        self.genome = None  # TODO: need way to determine genome before any posts occur!
        self.genome_warning = False
        self.annotation = None  # TODO: if appropriate, need way to determine annotation
        self.report_specs = {} # Unless the proper report specifications are discovered, a report can't be run
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
                if not self.genome_warning:
                    print >> sys.stderr, "WARNING: Mixing genomes in one report: %s and %s" % \
                                                        (self.genome, dxencode.GENOME_DEFAULTS[rep['organism']])
                    self.genome_warning  = True
                
        if verbose:
            print >> sys.stderr, "Replicates:"
            print >> sys.stderr, json.dumps(replicates,indent=4)
        return replicates


    def find_matching_files(self, files, out_type, rep_tech, sub_type=None, suffix=None):
        '''Returns list of matched files for out_type, replicate and possible ending.'''
        
        matching_files = []
        for file_obj in files:
            if file_obj.get('output_type') != out_type:
                continue
            if suffix != None:
                file_name = file_obj["submitted_file_name"]
                if not file_name.endswith(suffix):
                    continue
            if 'replicate' not in file_obj:
                print >> sys.stderr, "WARNING: %s %s has no 'replicate'" % \
                                    (file_obj['accession'],file_obj['submitted_file_name']) 
                continue
            br = file_obj['replicate']['biological_replicate_number']
            tr = file_obj['replicate']['technical_replicate_number']
            file_rep_tech = "rep%d_%d" % (br,tr)
            if file_rep_tech != rep_tech:
                continue
            acc = file_obj['accession'] 
            matching_files.append((out_type,sub_type,rep_tech,acc))
            self.obj_cache["exp"]["files"][acc] = file_obj
        return matching_files

    def get_relevant_files(self, exp_id, exp, report_specs, reps, verbose=False):
        '''Returns list of enc file_objects for selected files available on encoded.'''
        #verbose=True
        enc_files = []
        self.obj_cache["exp"]["files"] = {}
        
        files = dxencode.get_enc_exp_files(exp,report_specs["output_types"],lab="encode-processing-pipeline", \
                                                                                                key=self.server_key)
        # FIXME: special case to get around m2,m3
        new_list = []
        while len(files) > 0:
            file_obj = files.pop()
            if "genome_annotation" not in file_obj or file_obj["genome_annotation"] in [u'V19',u'M4']: # Note: this is a good 
               new_list.append(file_obj)
        files = new_list
        if len(files) == 0:
            if verbose:
                print >> sys.stderr, "Experiment "+exp_id+" has no results to report."
            return files

        # How to place them in order?
        for rep in reps:
            for out_type in report_specs["output_types"]:
                if "ending" in report_specs:
                    for sub_out, sub_type, suffix in report_specs["ending"]: 
                        if sub_out != out_type:
                            continue
                        matching_files = self.find_matching_files(files, out_type, rep['rep_tech'], sub_type, suffix)
                        if len(matching_files) > 0:
                            enc_files.extend(matching_files)
                else:
                    matching_files = self.find_matching_files(files, out_type, rep['rep_tech'])
                    if len(matching_files) > 0:
                        enc_files.extend(matching_files)
    
        if verbose:
            print >> sys.stderr, "Encoded files:"
            for (out_type,sub_type, rep_tech,acc) in enc_files:
                print >> sys.stderr, "'%s' %s %s %s" % \
                        (out_type,rep_tech,acc,self.obj_cache["exp"]["files"][acc]['submitted_file_name'])
        return enc_files


    def duration_from_enc_times(self,beg_time,end_time,second=False):
        '''Returns a duration string difference two encodeD formatted strings.'''
        if beg_time == None or end_time == None:
            return None
        beg_dt = datetime.strptime(beg_time, '%Y-%m-%dT%H:%M:%SZ')
        end_dt = datetime.strptime(end_time, '%Y-%m-%dT%H:%M:%SZ')
        duration = (end_dt - beg_dt)
        return duration.total_seconds()


    def file_size_string(self,size):
        '''Returns a file size string in G, M denomiations or else the number itself.'''
        if size >= 1000000000:
            return "%.2fG" % (size / 1000000000.0)
        if size >= 1000000:
            return "%.2fM" % (size / 1000000.0)
        if size >= 100000:
            return "%.1fK" % (size / 1000.0)
        return  size


    def get_special_metric(self,metric_id,file_obj,metric_key,metric_def,verbose=False):
        '''Returns a mocked up 'metric' object from specialized definitions.'''
        #verbose=True
        if metric_key != "file_cost":
            return None
        
        metric = {}
        notes = None
        file_notes = file_obj.get("notes")
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
                        step_notes = file_obj["step_run"].get("notes")
                        if step_notes != None:
                            try:
                                cost = json.loads(step_notes).get("dx_cost")
                            except:
                                cost = None
                    if cost != None:
                        metric[col] = cost
            elif col == "Time":
                #if notes != None:
                #    duration = notes.get("duration")
                step_run = self.enc_lookup_json(file_obj["step_run"])
                if step_run != None:
                    step_details = step_run["dx_applet_details"][0]
                    duration = self.duration_from_enc_times(step_details.get("started_running"),step_details.get("stopped_running"))
                    if duration != None:
                        metric[col] = duration
            elif col == "File Size":
                metric[col] = file_obj["file_size"]
        if verbose:
            print >> sys.stderr, "Special metric '"+metric_id+"':"
            print >> sys.stderr, json.dumps(metric,indent=4)
        return metric
            

    def get_metrics(self,exp_id,enc_files,report_specs,verbose=False):
        '''Returns a list of tuples of metrics objects for reporting statistics.'''
        #verbose=True
        # TODO: mock up a "metric that is not from a quality_metric object but from file and step_run info!
        metrics = []
        combo_metrics = []
        metric_ids = []
        self.obj_cache["exp"]["metrics"] = {}
        for (out_type,sub_type,rep_tech,acc) in enc_files:
            file_obj = self.obj_cache["exp"]["files"][acc]
            if sub_type != None:
                metric_keys = report_specs["metrics"].get(sub_type)
            else:
                metric_keys = report_specs["metrics"].get(out_type)
            
            assert metric_keys != None
            if not isinstance(metric_keys,list):
                metric_keys = [ metric_keys ]
            
            for metric_key in metric_keys:                
                metric_defs = self.METRIC_DEFS[metric_key]
                if "special_metric" in metric_defs and metric_defs["special_metric"]:
                    if verbose:
                        print >> sys.stderr, "Getting special metric for %s %s %s" % (exp_id,rep_tech,metric_key)
                    metric_id = acc + '/' + metric_key 
                    if metric_id not in metric_ids:  # Note that this will automatically give only one combined metric from 2 reps.
                        metric_ids.append(metric_id) 
                        metric = self.get_special_metric(metric_id,file_obj,metric_key,metric_defs)
                        if metric != None:
                            if metric_defs["per"] == "experiment":
                                combo_metrics.append((metric_key,"combined",metric_id))
                            else:
                                metrics.append((metric_key,rep_tech,metric_id))
                            self.obj_cache["exp"]["metrics"][metric_id] = metric
                else:
                    if verbose:
                        print >> sys.stderr, "Getting qc metric for %s %s %s" % (exp_id,rep_tech,metric_key)
                    file_qc_metrics = file_obj.get('quality_metrics')
                    if file_qc_metrics == None or len(file_qc_metrics) == 0:
                        if verbose:
                            print >> sys.stderr, "Experiment "+exp_id+" file "+acc+" has no quality_metric objects."
                        continue
                    for metric in file_qc_metrics:
                        #print >> sys.stderr, json.dumps(metric,indent=4)
                        #if verbose:
                        #    print >> sys.stderr, "  examining %s" % (metric["@type"])
                        if "name" in metric_defs:
                            if metric_defs["name"]+"_quality_metric" not in metric["@type"]:
                                continue
                        elif metric_key+"_quality_metric" not in metric["@type"]:
                            continue
                        metric_id = metric["@id"]
                        if metric_id not in metric_ids:  # Note that this will automatically give only one combined metric from 2 reps.
                            metric_ids.append(metric_id) 
                            if metric_defs["per"] == "experiment":
                                combo_metrics.append((metric_key,"combined",metric_id))
                            else:
                                metrics.append((metric_key,rep_tech,metric_id))
                            self.obj_cache["exp"]["metrics"][metric_id] = metric
        if len(combo_metrics) > 0:
            metrics.extend(combo_metrics)

        if verbose:
            print >> sys.stderr, "Encoded metrics:"
            for (metric_key,rep_tech,metric_id) in metrics:
                print >> sys.stderr, "%s %s %s" % (metric_key,rep_tech,metric_id)
        return metrics


    def type_headers(self,section_type,report_specs):
        '''Generates the headers for a type section from report specs.'''
        type_header = ""
        metric_keys = report_specs["metrics"][section_type]
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
            sub_type_found = False
            if "ending" in report_specs:
                for sub_out, sub_type, suffix in report_specs["ending"]:
                    if sub_out != out_type:
                        continue
                    Header += self.type_headers(sub_type,report_specs)
                    sub_type_found = True
            if not sub_type_found:
                Header += self.type_headers(out_type,report_specs)
            
        print Header

    def init_type_totals(self,section_type,report_specs):
        '''Intitializes totals for a type section from report specs.'''
        metric_keys = report_specs["metrics"][section_type]
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
        for out_type in report_specs["output_types"]:
            sub_type_found = False
            if "ending" in report_specs:
                for sub_out, sub_type, suffix in report_specs["ending"]:
                    if sub_out != out_type:
                        continue
                    self.init_type_totals(sub_type,report_specs)
                    sub_type_found = True
            if not sub_type_found:
                self.init_type_totals(out_type,report_specs)

    def format_value(self,val,format_type):
        '''Returns formatted string.'''
        if format_type == "file_size":
            return self.file_size_string(val)
        elif format_type == "duration":
            return dxencode.duration_string(val,include_seconds=False)
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
            
            
    def print_metrics(self,exp_id,enc_files,metrics,report_specs, verbose=False):
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

        # counts:
        print "# Experiments:\t" + str(self.exp_count)
        print "# Replicates:\t" + str(self.rep_count)


    def run(self):
        '''Runs mission_log from start to finish using command line arguments.'''
        args = self.get_args()
            
        self.server_key = args.server
        self.authid, self.authpw, self.server = dxencode.processkey(self.server_key)
        
        # Look up report type
        if args.report_type not in self.REPORTS_SUPPORTED or args.report_type not in self.REPORT_SPECS.keys():
            print >> sys.stderr, "Report type %s is not supported" % (args.report_type)
            sys.exit(1)

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
            self.exp = dxencode.get_exp(self.exp_id,must_find=True,key=self.server_key)
            if self.exp == None or self.exp["status"] == "error":
                print >> sys.stderr, "Unable to locate experiment %s in encoded (%s)" % (self.exp_id, self.server_key)
                continue
                
            # Type of experiment?
            exp_type = dxencode.get_assay_type(self.exp_id,self.exp)
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
            
            # Get files from encoded
            self.reps = self.find_replicates(self.exp_id, self.exp)
            enc_files = self.get_relevant_files(self.exp_id, self.exp, self.report_specs, self.reps,verbose=args.verbose)
            if enc_files == None or len(enc_files) == 0:
                print >> sys.stderr, "WARNING: No files available on encodeD for %s" % self.exp_id
                continue
            
            # Get metrics per file
            metrics = self.get_metrics(self.exp_id,enc_files,self.report_specs,verbose=args.verbose)
            if metrics == None or len(metrics) == 0:
                print >> sys.stderr, "WARNING: No metrics available on encodeD for %s" % self.exp_id
                continue
            
            # Print metrics
            self.print_metrics(self.exp_id,enc_files,metrics,self.report_specs)
            
        # print any totals available
        if self.report_specs:     
            self.print_totals()
        print >> sys.stderr, "(finished)"


if __name__ == '__main__':
    '''Run from the command line.'''
    mission_log = Mission_log()
    mission_log.run()

