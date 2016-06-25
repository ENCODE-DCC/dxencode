#!/usr/bin/env python2.7
# rollover.py 1.0.0
#
# Rollover.py will move all fastq files from hg19 to GRCh38 directiries

import argparse,os, sys
import json, urlparse, subprocess, itertools, logging, time
from datetime import datetime
from base64 import b64encode
import commands

import dxpy
import dx
import encd

class Rollover(object):
    '''
    Rollover module moves fastq files from hg19 to GRCh38 directories in dnanexus.
    '''
    TOOL_IS = 'rollover'
    HELP_BANNER = "Rollover moves fastq files from hg19 to GRCh38. " 
    ''' This help banner is displayed by get_args.'''

    SERVER_DEFAULT = 'www'
    '''This the server to makes files have been posed to.'''

    FOLDER_DEFAULT = "/"
    '''Where to start the search for experiment folders.'''

    EXPERIMENT_TYPES_SUPPORTED = [ 'long-rna-seq', 'small-rna-seq', 'rampage', 'dna-me' ] #, 'dnase' ]
    '''This module supports only these experiment (pipeline) types.'''
    
    ASSEMBLIES_SUPPORTED = { "hg19": "hg19", "GRCh38": "GRCh38", "mm10": "mm10" }
    '''This module supports only these assemblies.'''

    ANNOTATIONS_SUPPORTED = [ 'V24', 'V19', 'M2', 'M3', 'M4' ]
    '''This module supports only these annotations.'''
    
    def __init__(self):
        '''
        Rollover expects one or more experiment ids as arguments and will find fastq files that 
        should be moved from hg19 to GRCh38 directories.
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
        #self.annotation = None  # if appropriate (mice), points the way to the sub-dir
        #self.pipeline = None # pipeline definitions (filled in when experiment type is known)
        self.replicates = None # lost replicate folders currently found beneath experiment folder
        self.test = True # assume Test until told otherwise
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

        ap.add_argument('--verbose',
                        help='More debugging output.',
                        action='store_true',
                        required=False)

        if parse:
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
            print "Requested rolling %d experiments" % len(sorted_exp_ids)
            return sorted_exp_ids
            
        return []


    def find_fastq_files(self,exp_folder,replicates,verbose=False):
        '''Returns tuple list of (rep_tech,fid) of fastq files found in replicate folders.'''
        #verbose=True
        fastq_files = []
        for rep_tech in replicates:
            fids = []
            for glob in ['*.fastq.gz','*.fq.gz','*.fastq','*.fq']:
                some_fids = dx.find_file(exp_folder + rep_tech + '/' + glob,self.proj_id, multiple=True,recurse=False)
                if some_fids != None and len(some_fids) > 0:
                    if verbose:
                        print >> sys.stderr, " - Found %d matching '%s'" % (len(some_fids), glob)
                    fids.extend( some_fids )
            if verbose:
                print >> sys.stderr, " - Found %d fastq files" % (len(fids))
            if len(fids) > 0: 
                for fid in fids:
                    acc = dx.file_get_property("accession",fid)
                    if acc != None:
                        if verbose:
                            print >> sys.stderr, " - Verified by accession %s" % (acc)
                        fastq_files.append( (rep_tech,fid) )
                    else:
                        file_root = dx.description_from_fid(fid)['name'].split('.')[0]
                        if len(file_root) == 11:   # Expecting ENCFF000ABC.fastq.gz
                            if verbose:
                                print >> sys.stderr, " - Verified by file_root %s" % (file_root)
                        fastq_files.append( (rep_tech,fid) )

        if verbose:
            print >> sys.stderr, "Fastq files:"
            print >> sys.stderr, json.dumps(fastq_files,indent=4)
        return fastq_files

    def rep_fids_from_fastq_list(self,rep_tech,fastq_files):
        """Returns a simple list of fids that are associated with the replicate specified."""
        fids = []
        for (source_rep,fid) in fastq_files:
            if source_rep == rep_tech:
                fids.append(fid)
        return fids

    def run(self):
        '''Runs rollover from start to finish using command line arguments.'''
        args = self.get_args()
        self.test = args.test
            
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
        total_moved = 0
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
            self.hg19_umbrella_folder = dx.umbrella_folder(args.folder,self.FOLDER_DEFAULT,self.proj_name,self.exp_type,"posted","hg19")
            self.GRCh38_umbrella_folder = dx.umbrella_folder(args.folder,self.FOLDER_DEFAULT,self.proj_name,self.exp_type,"posted","GRCh38")
            if args.test:
                print "- Umbrella folders: hg19   '%s'" % (self.hg19_umbrella_folder)
                print "- Umbrella folders: GRCh38 '%s'" % (self.GRCh38_umbrella_folder)
            self.hg19_exp_folder = dx.find_exp_folder(self.project,exp_id,self.hg19_umbrella_folder,warn=True)
            if self.hg19_exp_folder == None:
                continue
            self.GRCh38_exp_folder = dx.find_exp_folder(self.project,exp_id,self.GRCh38_umbrella_folder,warn=True)
            if self.GRCh38_exp_folder == None:
                continue
            exp_count += 1
            print "- Potentially moving %s to '%s'..." % \
                                            (self.hg19_exp_folder, self.GRCh38_exp_folder)

            # 3) Given the experiment type, determine the expected results
            self.hg19_replicates = dx.find_replicate_folders(self.project,self.hg19_exp_folder, verbose=args.verbose)
            self.GRCh38_replicates = dx.find_replicate_folders(self.project,self.GRCh38_exp_folder, verbose=args.verbose)

            # 4) Get the list of fastq files in the hg19 directory
            fastq_files = self.find_fastq_files(self.hg19_exp_folder, self.hg19_replicates, verbose=args.verbose)
            print "- Found %d files that are available to move." % len(fastq_files)
            if len(fastq_files) == 0:
                continue

            # 5) For each hg19 replicate:
            files_moved = 0
            for rep_tech in self.hg19_replicates:
                sys.stdout.flush() # Slow running job should flush to piped log
                rep_fids = self.rep_fids_from_fastq_list(rep_tech,fastq_files)
                if len(rep_fids) > 0:
                    destination_folder = self.GRCh38_exp_folder + rep_tech + "/"
                    if self.test:
                        print "  * Would move %d file(s) to %s..." % (len(rep_fids),destination_folder)
                    else:
                        print "  * Moving %d file(s) to %s..." % (len(rep_fids),destination_folder)
                        dx.move_files(rep_fids,destination_folder,self.proj_id)
                    files_moved += len(rep_fids) 

            if not args.test:
                print "- For %s processed %d fasqs(s), moved %d file(s)" % (self.exp_id, len(fastq_files), files_moved)
            else:
                print "- For %s processed %d fastq(s), would move %d file(s)" % (self.exp_id, len(fastq_files), files_moved)
            total_moved += files_moved
            
        if not args.test:
            print "Processed %d experiment(s), moved %d file(s)" % \
                                                      (exp_count, total_moved)
        else:
            print "Processed %d experiment(s), would move %d file(s)" % \
                                                      (exp_count, total_moved)
        print "(finished)"


if __name__ == '__main__':
    '''Run from the command line.'''
    rollover = Rollover()
    rollover.run()

