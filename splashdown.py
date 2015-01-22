#!/usr/bin/env python
# splashdown.py 0.0.1
#
# Initial starting point accessonator.py in tf_chipseq.py and lrnaSplashdown.py 
#
# Splasdown is meant to run outside of dnanexus and to examine experiment directories to 
# find results up upload to encoded.
#
# 1) Lookup experiment type from encoded, based on accession
# 2) Locate the experiment accession named folder
# 3) Given the experiment type, determine the expected results
# 4) Given expected results locate any files (by glob) that should be uploaded for 
#    a) each single replicate (in replicate sub-folders named as reN_N/
#    b) combined replicates in the experiment folder itself
# 5) For each file that should be uploaded, determine if the file needs to be uploaded (not already uploaded).
# 6) For each file that needs to be uploaded:
#    a) discover all necessary dx information needed for upload.
#    b) gather any other information necessary from dx and encoded. (notice hand-waving)
#    c) Upload file and update encoded database. 
#    d) Update dnanexus file with file accession tag.
# 7) Either exit or advance to the next experiment folder

import argparse,os, sys
import json, urlparse, subprocess, itertools
#import requests, logging, re, shlex, time
#from datetime import datetime

import dxpy
import dxencode


class Splashdown(object):
    PROJECT_DEFAULT = 'scratchPad'
    '''This the default DNA Nexus project to use for the long RNA-seq pipeline.'''
    
    SERVER_DEFAULT = 'test'
    '''This the default server to post files to.'''
    
    RESULT_FOLDER_DEFAULT = "/"
    '''Where to start the search for experiment folders.'''

    ASSEMBLIES_SUPPORTED = { "hg19": "GRCh37", "hg38": "GRCh37", "mm10": "GRCm38" }
    '''This module supports only these assemblies.'''

    ANNOTATIONS_SUPPORTED = [ 'v19', 'M2', 'M3', 'M4' ]
    '''This module supports only these annotations.'''
    
    EXPERIMENT_TYPES_SUPPORTED = [ "long-rna-seq", "small-rna-seq", "rampage" ] #,"dna-me","chip-seq" ]
    '''This module supports only these experiment (pipeline) types.'''

    #Pipeline specifications include order of steps, steps per replicate, combined steps and 
    #within steps, the output_type: file_glob that define expected results.
    PIPELINE_SPECS = {
        "long-rna-seq": {
            "step-order": [ "align-tophat","signals-top-se","signals-top-pe",
                            "align-star","signals-star-se","signals-star-pe","quant-rsem"],
            "combined":   {},
            "replicate":  {
                "align-tophat":    { "tophat_bam":           "*_tophat.bam"                 },
                "signals-top-se":  { "tophat_all_bw":        "*_tophat_all.bw",
                                     "tophat_uniq_bw":       "*_tophat_uniq.bw"             },
                "signals-top-pe":  { "tophat_minus_all_bw":  "*_tophat_minusAll.bw",
                                     "tophat_minus_uniq_bw": "*_tophat_minusUniq.bw",
                                     "tophat_plus_all_bw":   "*_tophat_plusAll.bw",
                                     "tophat_plus_uniq_bw":  "*_tophat_plusUniq.bw"         },
                "align-star":      { "star_genome_bam":      "*_star_genome.bam",
                                     "star_anno_bam":        "*_star_anno.bam"              },
                "signals-star-se": { "star_all_bw":          "*_star_genome_all.bw",
                                     "star_uniq_bw":         "*_star_genome_uniq.bw"        },
                "signals-star-pe": { "star_minus_all_bw":    "*_star_genome_minusAll.bw",
                                     "star_minus_uniq_bw":   "*_star_genome_minusUniq.bw",
                                     "star_plus_all_bw":     "*_star_genome_plusAll.bw",
                                     "star_plus_uniq_bw":    "*_star_genome_plusUniq.bw"    },
                "quant-rsem":      { "rsem_iso_results":     "*_rsem.isoforms.results",
                                     "rsem_gene_results":    "*_rsem.genes.results"         }
            }
        },
        "small-rna-seq": {
            "step-order": [ "align","signals"],
            "combined":   {},
            "replicate":  {
                "align":           { "genome_bam":           "*_star_genome.bam"            },
                "signals":         { "all_plus_bw":          "*_small_plusAll.bw",
                                     "all_minus_bw":         "*_small_minusAll.bw",
                                     "unique_plus_bw":       "*_small_plusUniq.bw",
                                     "unique_minus_bw":      "*_small_minusUniq.bw"         }
            }
        },
        "rampage": {
            "step-order": [ "align","signals","peaks","idr"],
            "combined":   {
                "idr":             { "rampage_idr_gff":      "*_rampage_idr.gff",
                                     "rampage_idr_bb":       "*_rampage_idr.bb" }
            },
            "replicate":  {
                "align":           { "rampage_marked_bam":   "*_rampage_star_marked.bam" },
                "signals":         { "all_plus_bw":          "*_rampage_5p_plusAll.bw",    
                                     "all_minus_bw":         "*_rampage_5p_minusAll.bw",
                                     "unique_plus_bw":       "*_rampage_5p_plusUniq.bw", 
                                     "unique_minus_bw":      "*_rampage_5p_minusUniq.bw" },
                "peaks":           { "rampage_peaks_bb":     "*_rampage_peaks.bb",
                                     "rampage_peaks_gff":    "*_rampage_peaks.gff" }
            }
        }
    }

    def __init__(self):
        '''
        Splashdown expects one or more experiment ids as arguments and will find, document
        and upload files in the associated directory. 
        '''
        self.args = {} # run time arguments
        self.server_key = 'test'
        self.proj_name = self.PROJECT_DEFAULT
        self.project = None
        self.proj_id = None
        self.exp = {}  # Will hold the encoded exp json
        self.exp_id = None
        self.exp_type = {}  # Will hold the experiment's assay_type, normalized to known tokens.
        self.genome = None  # TODO: need way to determine genome before any uploads occur!
        self.annotation = None  # TODO: if appropriate, need way to determine annotation
        self.pipeline = None # pipeline definitions (filled in when experiment type is known)
        self.replicates = None # lost replicate folders currently found beneath experiment folder
        self.test = True # assume Test until told otherwise
    
    def get_args(self):
        '''Parse the input arguments.'''
        ### PIPELINE SPECIFIC
        ap = argparse.ArgumentParser(description="Handles splashdown of launched pipeline runs " +
                    "for supported experiment types. Can be run repeatedly and will only try to " +
                    "upload result files that have not been previously uploaded. All results " +
                    "are expected to be in folder /<resultsLoc>/<experiment> and any replicate " +
                    "sub-folders named as " +
                    "<experiment>/rep<biological-replicate>_<technical-replicate>.")
        ### PIPELINE SPECIFIC

        ap.add_argument('-e', '--experiments',
                        help='One or more ENCODED experiment accessions',
                        nargs='+',
                        required=True)

        #ap.add_argument('--br', '--biorep',
        #                help="Biological Replicate number (default: 1)",
        #                type=int,
        #                default='1',
        #                required=False)

        #ap.add_argument('--tr', '--techrep',
        #                help="Technical replicate number (default: 1)",
        #                type=int,
        #                default='1',
        #                required=False)

        #ap.add_argument('--cr','--combine-replicates',
        #                help="Combine or compare two replicates (e.g.'1 2_2').'",
        #                nargs='+',
        #                required=False)

        ap.add_argument('--project',
                        help="Project to run analysis in (default: '" + self.PROJECT_DEFAULT + "')",
                        default=self.PROJECT_DEFAULT,
                        required=False)

        #ap.add_argument('--refLoc',
        #                help="The location to find reference files (default: '" + \
        #                                    REF_PROJECT_DEFAULT + ":" + REF_FOLDER_DEFAULT + "')",
        #                default=REF_FOLDER_DEFAULT,
        #                required=False)

        ap.add_argument('--results_folder',
                        help="The location to to place results folders (default: '<project>:" + \
                                                                self.RESULT_FOLDER_DEFAULT + "')",
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

    
    def find_exp_folder(self,exp_id,results_folder='/'):
        '''Finds the experiment folder, given an accession.'''
        # normalize
        if not results_folder.startswith('/'):
            results_folder = '/' + results_folder
        if not results_folder.endswith('/'):
            results_folder += '/'
        target_folder = dxencode.find_folder(exp_id,self.project,results_folder)
        if target_folder == None or target_folder == "":
            print "Unable to locate target folder for %s in project %s" % (exp_id, self.proj_name)
            return None
        #return self.proj_name + ':' + target_folder + '/'
        return target_folder + '/'
        

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
        self.genome = None  # TODO: need way to determine genome before any uploads occur!
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
        if ext in ["bam","bw","bb","gff"]:
            return ext
        if ext == "results":
            return "csv"
        elif ext == "bw":
            return "bigwig"
        elif ext == "bb":
            return "bigbed"
        return None


    def find_genome(self,path):
        '''Try to determine genome from file path.'''
        # TODO: There should be a better way of determining genome and annotation!!!
        if self.genome == None:
            for part in path.split('/'):
                if part in self.ASSEMBLIES_SUPPORTED.keys():
                    self.genome = self.ASSEMBLIES_SUPPORTED[part]
                elif self.annotation == None and part in self.ANNOTATIONS_SUPPORTED:
                    self.annotation = part
        print "- Genome: %s, annotation: %s" % (self.genome,self.annotation)
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
                if self.genome == None:
                    self.genome = self.find_genome(dxencode.file_path_from_fid(fid))
            else:
                return []      # Only include files from completed steps!
        return step_files

    def find_expected_files(self,exp_folder,replicates,verbose=False):
        '''Returns tuple list of (type,rep_tech,fid) of files expected to be uploaded to ENCODE.'''
        expected = []
        # First find replicate step files
        for rep_tech in replicates:
            for step in self.pipeline["step-order"]:
                if step not in self.pipeline["replicate"]:
                    continue
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
        file_name = dxencode.file_path_from_fid(fid,projectToo=True)
        file_size = dxencode.description_from_fid(fid).get('size')

        # TODO: encoded should allow searching on fid !!!
        authid, authpw, server = dxencode.processkey(self.server_key)
        url = urlparse.urljoin(server,
                'search/?type=file&submitted_file_name=%s&format=json&frame=object' % file_name)
        response = dxencode.encoded_get(url,authid,authpw)
        try:
            if verbose:
                print "-- Looking for %s (size:%d) in %s" % (file_name, file_size, server)
            response.raise_for_status()
            if response.json()['@graph']:
                encode_file = response.json()['@graph'][0]
                #logger.info("Found potential duplicate: %s" %(duplicate_item.get('accession')))
                if file_size ==  encode_file.get('file_size'):
                    #logger.info("%s %s: File sizes match, assuming duplicate." % \
                    #                               (str(file_size), encode_file.get('file_size')))
                    return encode_file
                else:
                    #logger.info("%s %s: File sizes differ, assuming new file." % \
                    #                                (str(file_size), encode_file.get('file_size')))
                    return None
            else:
                #logger.info("No duplicate ... proceeding")
                return None
        except:
            #logger.warning('Duplicate accession check failed: %s %s' % (response.status_code, response.reason))
            #logger.debug(response.text)
            return None

    def find_needed_files(self,files_expected,verbose=False):
        '''Returns the tuple list of files that NEED to be uploaded to ENCODE.'''
        needed = []
        for (out_type, rep_tech, fid) in files_expected:
            # Don't bother searching for tag: tag means upload requested, NOT upload succeeded!
            #for tag in dxencode.description_from_fid(fid)['tags']:
            #    m = re.search(r'(ENCFF\d{3}\D{3})|(TSTFF\D{6})', tag)

            if self.find_in_encode(fid,verbose) == None:
                needed.append( (out_type,rep_tech,fid) )
        if verbose:
            print "Needed files:"
            print json.dumps(needed,indent=4)
        return needed


    def find_derived_from(self,job,verbose=False):
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
            inp_obj = dxencode.description_from_fid(inp_fid)
            if inp_obj != None:
                if self.genome == None:
                    self.genome = self.find_genome(inp_obj["folder"])
                if "properties" in inp_obj and "accession" in inp_obj["properties"]:
                    derived_from.append(inp_obj["properties"]["accession"])
                else: # if file name is primary input (fastq) and is named as an accession
                    if inp_obj["name"].startswith("ENCFF"):
                        if inp_obj["name"].endswith(".fastq.gz") \
                        or inp_obj["name"].endswith(".fq.gz") \
                        or inp_obj["name"].endswith(".fastq") \
                        or inp_obj["name"].endswith(".fq"):
                            root = inp_obj["name"].split('.')[0]
                            for acc in root.split('-'): #usually only one
                                # TODO: make sure concat files do what they are supposed to!!!
                                if acc.startswith("ENCFF") and len(acc) == 11:
                                    derived_from.append(acc)
        if verbose:
            print "Derived files:"
            print json.dumps(derived_from,indent=4)
        return derived_from


    def add_encoded_info(self,obj,rep_tech,fid,verbose=False):
        '''Updates an object with information from encoded database.'''
        obj['lab'] = self.exp['lab']['@id']
        obj['award'] = self.exp['award']['@id']
        
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


    def make_upload_obj(self,out_type,rep_tech,fid,verbose=False):
        '''Returns an object for submitting a file to encode, with all dx info filled in.'''
        upload_obj = {}
        upload_obj['dataset'] = self.exp_id
        upload_obj["output_type"] = out_type
        dx_obj = dxencode.description_from_fid(fid)
        #if verbose:
        #    print "dx_obj:"
        #    print json.dumps(dx_obj,indent=4)
        job = dxencode.job_from_fid(fid)
        #if verbose:
        #    print "job:"
        #    print json.dumps(job,indent=4)
        #applet = dxencode.applet_from_fid(fid)
        upload_obj["file_format"] = self.file_format(dx_obj["name"])
        upload_obj["derived_from"] = self.find_derived_from(job)
        upload_obj['submitted_file_name'] = dxencode.file_path_from_fid(fid,projectToo=True)
        upload_obj['file_size'] = dx_obj["size"]
        #upload_obj['md5sum'] = calculated_md5 # TODO: Find from file properties???
        if self.genome == None:
            print "Error: could not determine genome assembly!"
            sys.exit(1)
        upload_obj['assembly'] = self.genome
        if self.annotation != None:
            upload_obj['genome_annotation'] = self.annotation

        dxFile = dxencode.file_handler_from_fid(fid)
        versions = dxencode.get_sw_from_log(dxFile, '\* (\S+)\s+version:\s+(\S+)')
        notes = dxencode.create_notes(dxFile, versions)
        #notes['qc'] = flagstat_parse(bamqc) ????
        # TODO: find json added to job as a result of returns? ??
        if "notes" in job["output"]:
           notes.update(json.load(job["output"]["json"]))
        upload_obj['notes'] = json.dumps(notes)

        #print "  - Adding encoded information."
        upload_obj = self.add_encoded_info(upload_obj,rep_tech,fid) 

        if verbose:
            print "Upload_obj from dx info:"
            print json.dumps(upload_obj,indent=4)
        return upload_obj


    def file_upload(self,upload_obj,test=True):
        '''Uploads a file to encoded.'''
        path = upload_obj['submitted_file_name'].split(':')[1]
        if test:
            print "  - Test upload %s to '%s' server" % (path,self.server_key)
            return "ENCFF00FAKE"
        else:
            (AUTHID,AUTHPW,SERVER) = processkey(self.server_key)
            try:
                result = dxencode.encoded_post_file(upload_obj,SERVER,AUTHID,AUTHPW)
                print "  - Real upload %s to '%s' server" % (path,self.server_key)
                return result.get('accession')
            except:
                print "  * FAILED upload of %s to '%s' server" % (path,self.server_key)
        return None


    def file_mark_accession(self,fid,accession,test=True):
        '''Adds/replaces accession to a file's properties.'''
        file_handler = dxencode.file_handler_from_fid(fid)
        properties = file_handler.get_properties()
        path = dxencode.file_path_from_fid(fid)
        if "accession" in properties and properties["accession"] != accession:
            print "Warning: file %s has accession %s but is being uploaded as accession %s" % \
                (path,properties["accession"],accession)
            #sys.exit(1)
        properties["accession"] = accession
        if test:
            print "  - Test update %s with accession='%s'" % (path,accession)
        else:
            file_handler.set_properties(properties)
            print "  - Real update %s with accession='%s'" % (path,accession)


    def run(self):
        '''Runs splasdown from start to finish using command line arguments.'''
        args = self.get_args()
        self.server_key = args.server
        self.proj_name = args.project
        self.project = dxencode.get_project(self.proj_name)
        self.proj_id = self.project.get_id()
        print "== Running in project [%s] and will post to the [%s] server ==" % \
                                                        (self.proj_name,self.server_key)
        
        for exp_id in args.experiments:
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
            self.exp_folder = self.find_exp_folder(exp_id,args.results_folder)
            if self.exp_folder == None:
                print "Unable to locate experiment folder %s in dnanexus" % exp_id
                continue
            print "- Examining %s:%s for '%s' results..." % \
                                            (self.proj_name, self.exp_folder, self.exp_type)

            # 3) Given the experiment type, determine the expected results
            self.pipeline   = self.pipeline_specification(args,self.exp_type,self.exp_folder)
            self.replicates = self.find_replicate_folders(self.exp_folder)

            # 4) Given expected results locate any files (by glob) that should be uploaded for 
            #    a) each single replicate (in replicate sub-folders named as reN_N/
            #    b) combined replicates in the experiment folder itself
            files_expected = self.find_expected_files(self.exp_folder,self.replicates)
            print "- Found %d files that are available to upload." % len(files_expected) 
            if len(files_expected) == 0:
                continue

            # 5) For each file that should be uploaded, determine if the file needs to be uploaded (not already uploaded).
            files_to_upload = self.find_needed_files(files_expected)
            print "- Found %d files that need to be uploaded" % len(files_to_upload) 
            if len(files_to_upload) == 0:
                continue

            # 6) For each file that needs to be uploaded:
            for (out_type,rep_tech,fid) in files_to_upload:
                # a) discover all necessary dx information needed for upload.
                # b) gather any other information necessary from dx and encoded. (notice hand-waving)
                print "  Document file %s" % dxencode.file_path_from_fid(fid) 
                upload_obj = self.make_upload_obj(out_type,rep_tech,fid)

                # c) Upload file and update encoded database. 
                accession = self.file_upload(upload_obj,args.test)

                # d) Update dnanexus file with file accession tag.
                if accession != None:
                    self.file_mark_accession(fid,accession,args.test)

        print "(finished)"
                

if __name__ == '__main__':
    '''Run from the coomad line.'''
    recovery = Splashdown()
    recovery.run()

