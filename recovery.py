#!/usr/bin/env python2.7
# recovery.py 0.0.1
#
# Recovery is derived from splashdow.py but is meant to update (older) already posted files to 
# fill in missing info not handled properly when the files were originally posted.
#
# In particular:
# 1) Fill in in 'analysis_step_runs' and 'workflow_runs', where missing.sys.stderr.write(
# 2) Fill in derived_from for reference files not included in the original posting.
# 3) TODO: eventually fill in QC_metrics (though that may be handled by splashdown itself.

import argparse,os, sys
import json, urlparse, subprocess, itertools, logging, time
#import requests, re, shlex, time
from datetime import datetime

import dxpy
import dxencode
from splashdown import Splashdown

class Recovery(Splashdown):
    '''
    Recovery module updates posted files with information not handled during original post.
    Descendent from 'Splashdown' class.
    '''

    HELP_BANNER = "Handles recovery of missing information for already posted files. For all files already posted " + \
                  "for an experiment, will generate normal payload package from DX and attempt to update ENCODEd " + \
                  "with missing information. " 
    ''' This help banner is displayed by get_args.'''

    EXPERIMENT_TYPES_SUPPORTED = [ 'long-rna-seq' ] #, 'small-rna-seq', 'rampage' ] #,"dna-me","chip-seq" ]
    '''This module supports only these experiment (pipeline) types.'''

    def __init__(self):
        Splashdown.__init__(self)
        self.step_runs_patched = 0
        self.exp_files = None
        
    def get_args(self):
        '''Parse the input arguments.'''
        ap = Splashdown.get_args(self,parse=False)
        
        ap.add_argument('--start_at',
                        help="Start processing with this file accession.",
                        default=None,
                        required=False)
        return ap.parse_args()

    def find_in_encode_old_way(self,fid,exp_files,verbose=False):
        '''Looks for file in encode with the same submitted file name as the fid has.'''
        file_path = dxencode.file_path_from_fid(fid,projectToo=True)
        file_size = dxencode.description_from_fid(fid).get('size')
        file_name = file_path.split('/')[-1]
        
        if verbose:
            sys.stderr.write("Looking for '%s' of size: %d\n" % (file_name, file_size))

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
                sys.stderr.write("  %s %d\n" % (enc_file.get("submitted_file_name"),exp_file.get('file_size')))
                
                    
        if found_file != None and verbose:
            sys.stderr.write("Found file:\n")
            sys.stderr.write(json.dumps(found_file,indent=4) + '\n')
        return found_file
        return None

    def find_posted_files(self,files_expected,test=True,verbose=False):
        '''Returns the tuple list of files already posted to ENCODE.'''
        
        # get all files associated with the experiment up front, just in case it is needed:
        self.exp_files = dxencode.get_enc_exp_files(self.exp,key=self.server_key)
        if len(self.exp_files) == 0:
            print "* ERROR: found no files associated with experiment: " + self.exp_id

        posted = []
        self.found = {}
        for (out_type, rep_tech, fid) in files_expected:
            if fid in self.found.keys():
                continue  # No need to handle the same file twice
            if verbose:
                sys.stderr.write("* DEBUG Working on: " + dxencode.file_path_from_fid(fid,projectToo=True) + '\n')
                
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
                    if verbose:
                        sys.stderr.write("* DEBUG   Accession: " + accession + '\n')
                    f_obj =  self.enc_file_find_find_by_dxid(fid)  # Look by alias first incase ENCFF v. TSTFF mismatch
                    if f_obj != None: # Verifyably posted
                        posted.append( (out_type,rep_tech,fid) )
                        self.found[fid] = f_obj
                        # TODO: Compare accessions.  But not ENCFF to TSTFF
                        if verbose:
                            sys.stderr.write("* DEBUG   Found by alias\n")
                        continue

                    # look by accession
                    if verbose:
                        sys.stderr.write("* DEBUG   Not found by alias: dnanexus:" + fid + '\n')
                    f_obj = self.enc_lookup_json( 'files/' + accession,must_find=False)
                    if f_obj != None: # Verifyably posted
                        posted.append( (out_type,rep_tech,fid) )
                        self.found[fid] = f_obj
                        # update encoded file with alias
                        update_payload = {}
                        update_payload['aliases'] = [ 'dnanexus:' + fid ]
                        if not test:
                            ret = dxencode.encoded_patch_obj(accession, update_payload, self.server, self.authid, self.authpw)
                            #if ret == accession:
                            print "  * Updated ENCODEd '"+accession+"' with alias."
                        if verbose:
                            sys.stderr.write("* DEBUG   Found by accession: " + accession + '\n')
                        continue
                    else:
                        if verbose:
                            sys.stderr.write("* DEBUG   Not found by accession: " + accession + '\n')
                        
            if accession == '':
                if verbose:
                    print "* DEBUG Accession not found."
                # No accession in properties, but try to match by fid anyway.
                #f_obj =  self.find_in_encode(fid,verbose)
                f_obj =  self.enc_file_find_find_by_dxid(fid)
                if f_obj != None:
                    posted.append( (out_type,rep_tech,fid) )
                    self.found[fid] = f_obj
                    # update dx file property
                    accession = f_obj['accession']
                    if accession.startswith('ENCFF'):  # Only bother if it is the real thing.  Ambiguity with 'TSTFF'
                        ret = dxencode.dx_file_set_property(fid,'accession',accession,add_only=True,test=test)
                        if not test:
                            if ret == accession:
                                print "  * Updated DX property with accession."
                    if verbose:
                        sys.stderr.write("* DEBUG   Found by alias\n")
                    continue

                if verbose:
                    sys.stderr.write("* DEBUG   Not found by alias: dnanexus:" + fid + '\n')

            # Check the old fashion method of file name and size.
            f_obj =  self.find_in_encode_old_way(fid,self.exp_files) 
            if f_obj != None:
                posted.append( (out_type,rep_tech,fid) )
                self.found[fid] = f_obj
                # Update ENC and DX:
                accession = f_obj['accession']
                if accession.startswith('ENCFF'):  # Only bother if it is the real thing.  Ambiguity with 'TSTFF'
                    ret = dxencode.dx_file_set_property(fid,'accession',accession,add_only=True,test=test)
                    if not test:
                        if ret == accession:
                            print "  - Updated DX property with accession."
                        ret = dxencode.encoded_patch_obj(accession, update_payload, self.server, self.authid, self.authpw)
                        #if ret == accession:
                        print "  * Updated ENCODEd '"+accession+"' with alias."
                            
                if verbose:
                    sys.stderr.write("* DEBUG   Found the old fashion way.\n")
                continue
            if verbose:
                sys.stderr.write("* DEBUG   Not Found the old fashion way.\nVERBOSE...\n")
                f_obj =  self.find_in_encode_old_way(fid,self.exp_files,verbose=True) 

        if verbose:
            sys.stderr.write("Posted files:" + '\n')
            sys.stderr.write(json.dumps(posted,indent=4) + '\n')
        return posted

    def find_derived_from(self,fid,job,verbose=False):
        '''Returns list of accessions a file is drived from based upon job inouts.'''
        input_accessions = []
        input_file_count = 0
        #verbose=True
        for inp in job["input"].values():
            if not type(inp) == dict:
                continue # not a file input
            dxlink = inp.get("$dnanexus_link")
            if dxlink == None:
                continue
            if not type(dxlink) == dict:
                inp_fid = dxlink
            else:
                inp_fid = dxlink.get("id")
            inp_obj = dxencode.description_from_fid(inp_fid,properties=True)
            input_file_count += 1
            if verbose: # NOTE: verbose can help find the missing accessions
                sys.stderr.write("* derived from: " + inp_fid + " " + inp_obj["project"] + ":" + \
                                                                        dxencode.file_path_from_fid(inp_fid) + '\n')
            if inp_obj != None:
                self.genome = self.find_genome_annotation(inp_obj)
                acc_key = dxencode.dx_property_accesion_key('https://www.encodeproject.org') # Prefer production accession
                accession = None
                if "properties" in inp_obj:
                    #if acc_key not in inp_obj["properties"] and self.server_key != 'www':   # test, beta replicated from www
                    #    acc_key = dxencode.dx_property_accesion_key(self.server)
                    #    #acc_key = dxencode.dx_property_accesion_key('https://www.encodeproject.org')
                    if verbose:
                        sys.stderr.write("Input file properties... looking for '"+acc_key+"'\n") 
                        sys.stderr.write(json.dumps(inp_obj["properties"],indent=4) + '\n')
                    # older runs don't have accession in properties
                    if acc_key in inp_obj["properties"]:
                        accession = inp_obj["properties"][acc_key]
                        # Must check if file exists!!
                        file_obj = self.enc_lookup_json( 'files/' + accession,must_find=False)
                        if file_obj == None:
                            accession = None
                        else:
                            input_accessions.append(accession)
                    if accession == None and self.server_key != 'www':
                        acc_key = dxencode.dx_property_accesion_key(self.server)  # demote to server's accession
                        if verbose:
                            sys.stderr.write("Now looking for '"+acc_key+"'\n") 
                        if acc_key in inp_obj["properties"]:
                            accession = inp_obj["properties"][acc_key]
                            # Must check if file exists!!
                            file_obj = self.enc_lookup_json( 'files/' + accession,must_find=False)
                            if file_obj == None:
                                accession = None
                            else:
                                input_accessions.append(accession)
                if accession == None:
                    if verbose:
                        sys.stderr.write("Accession not found in properties or not in ENCODEd.\n") 
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
                if accession == None:
                    # If all else fails, try the old fashion way
                    if verbose:
                        sys.stderr.write("Accession still not found.  Trying the old fashion way.\n") 
                    f_obj =  self.find_in_encode_old_way(inp_fid,self.exp_files)
                    if f_obj != None and "accession" in f_obj:
                        accession = f_obj["accession"]
                        input_accessions.append(accession)

        if len(input_accessions) < input_file_count:
            if self.test:
                print "WARNING: not all input files are accounted for in 'derived_from'!"
                print json.dumps(input_accessions,indent=4)
            else:
                print "ERROR: not all input files are accounted for in 'derived_from'!"
                print json.dumps(input_accessions,indent=4)
                sys.exit(1)  # TODO: determine if this is appropriate
        
        # Now that we have the full derived_from, we can remove some         
        # UGLY special cases:
        derived_from = []
        for acc in input_accessions:
            if not acc.endswith("chrom.sizes"):  # Chrom.sizes is just a bit too much info
                derived_from.append(acc)
        
        if verbose:
            sys.stderr.write("Derived files: for " + dxencode.file_path_from_fid(fid) + '\n')
            sys.stderr.write(json.dumps(derived_from,indent=4) + '\n')
        return derived_from


    def pipeline_qualifiers(self,rep_tech,app_name=None):
        '''Determines and pipeline_qualifiers in ugly special-case code.'''
        #if app_name != None:
        #    print "Looking for pipe_qualifiers for '"+rep_tech+"' and '"+app_name+"'" # FIXME: debug
        #else:
        #    print "Looking for pipe_qualifiers for '"+rep_tech+"'" # FIXME: debug
        if "exp" in self.obj_cache and rep_tech in self.obj_cache["exp"] \
        and "pipe_qualifiers" in self.obj_cache["exp"][rep_tech]:
            #print "FOUND cached pipeline_qualifiers" # FIXME: debug
            return self.obj_cache["exp"][rep_tech]["pipe_qualifiers"]
        #if "exp" in self.obj_cache and pipe_qualifiers in self.obj_cache["exp"]:
        #    return self.obj_cache["exp"]["pipe_qualifiers"]
            
        pipe_qualifiers = {}
        pipe_qualifiers["version"] = "1"
        #if self.exp_type.startswith('long-rna-seq'): # FIXME: ugly special case
        #    pipe_qualifiers["version"] = "2"

        pipe_qualifiers["qualifier"] = ''
        if self.exp_type.startswith('long-rna-seq'): # FIXME: ugly special case
            if app_name and len(app_name):
                if app_name.endswith('-pe'):
                    pipe_qualifiers["qualifier"] = '-pe'
                elif app_name.endswith('-se'):
                    pipe_qualifiers["qualifier"] = '-se'
            if len(pipe_qualifiers["qualifier"]) == 0 and rep_tech and len(rep_tech):
                if rep_tech.startswith("rep") and not rep_tech.startswith("reps"):
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
        elif self.exp_type.startswith('small-rna-seq'): # FIXME: ugly special case
            pipe_qualifiers["qualifier"] = '-se'
        elif self.exp_type.startswith('rampage'): # FIXME: ugly special case
            pipe_qualifiers["qualifier"] = '-pe'
                            
        if len(pipe_qualifiers["qualifier"]) == 0:
            print "Error: Can't determine pipeline qualifiers for '%s' replicate of '%s' experiment" % \
                                                                                (replicate,self.exp_type)
        if "exp" not in self.obj_cache:
             self.obj_cache["exp"] = {}
        if rep_tech not in self.obj_cache["exp"]:
            self.obj_cache["exp"][rep_tech] = {}
        self.obj_cache["exp"][rep_tech]["pipe_qualifiers"] = pipe_qualifiers
        #self.obj_cache["exp"]["pipe_qualifiers"] = pipe_qualifiers
        return pipe_qualifiers


    def enc_step_run_find_or_create(self,job,dxFile,rep_tech,test=False,verbose=False):
        '''Finds or creates the 'analysis_step_run' encoded object that actually created the file.'''
        step_run = Splashdown.enc_step_run_find_or_create(self,job,dxFile,rep_tech,test=test,verbose=verbose)
        
        # If step_run was found in encoded, see if it needs to be updated:
        if '@id' not in step_run: # Must have just been created so no update required!
            return step_run
            
        # Do we need to patch step_run?
        update_step_run = {}
        patch_required = False
        
        # Alias   
        job_id = job.get('id')
        step_alias = 'dnanexus:' + job_id
        if 'aliases' not in step_run or step_alias not in step_run['aliases']:
            if 'aliases' not in step_run:
                step_run['aliases'] = []
            print "  + Need to update step_run['aliases'] new: '"+step_alias+"'"
            step_run['aliases'].append(step_alias)
            update_step_run['aliases'] = step_run['aliases']
            patch_required = True

        # Find workflow. Should be in cache already!
        # NOTE: not interested in updating wf_runs
        dx_wfr_id = job.get('analysis')
        wf_run = self.enc_workflow_run_find_or_create(dx_wfr_id,rep_tech,test=self.test,verbose=verbose)
        wf_run_notes = json.loads(wf_run["notes"])
        # The step_run["workflow_run"] will be object coming from enc.
        if "workflow_run" not in step_run:
            print "  + Need to update missing step_run['workflow_run'] new: '"+dx_wfr_id+"'"
            step_run["workflow_run"] = "/workflow-runs/dnanexus:" + dx_wfr_id
            update_step_run["workflow_run"] = "/workflow-runs/dnanexus:" + dx_wfr_id
            patch_required = True
        else:
            enc_wf = step_run["workflow_run"]
            if 'aliases' not in enc_wf:
                print "  + Need to update missing alias step_run['workflow_run'] new: '"+dx_wfr_id+"'"
                #step_run["workflow_run"] = "/workflow-runs/dnanexus:" + dx_wfr_id
                update_step_run["workflow_run"] = "/workflow-runs/dnanexus:" + dx_wfr_id
                patch_required = True
            else:
                patch_wf = True 
                for alias in enc_wf['aliases']:
                    if alias == "dnanexus:" + dx_wfr_id:
                        patch_wf = False
                        break
                if patch_wf:
                    print "  + Need to update step_run['workflow_run'] new: '"+dx_wfr_id+"'"
                    #print json.dumps(enc_wf['aliases'],indent=4,sort_keys=True)
                    #print json.dumps(enc_wf,indent=4,sort_keys=True)
                    #enc_wf['aliases'].append("/workflow-runs/dnanexus:" + dx_wfr_id)
                    #step_run["workflow_run"] = "/workflow-runs/dnanexus:" + dx_wfr_id
                    update_step_run["workflow_run"] = "/workflow-runs/dnanexus:" + dx_wfr_id
                    patch_required = True
            
        # Find analysis_step
        dx_app_id = job.get('applet')
        dx_app_name = job.get('executableName')
        dx_app_ver = self.find_app_version(dxFile)
        if dx_app_ver and 'version' in dx_app_ver:
            dx_app_ver = str( dx_app_ver.get('version') )
            if dx_app_ver[0] == 'v':
                dx_app_ver = dx_app_ver[1:]
        # FIXME: UGLY temporary special case!!!
        if dx_app_name  == "rampage-peaks" and dx_app_ver == "1.0.1":
            dx_app_version = "1.1.0"  # Because the version was supposed to bump the second digit if a tool changes.
        # FIXME: UGLY temporary special case!!!
        if not dx_app_ver or not isinstance(dx_app_ver, str) or len(dx_app_ver) == 0:
            print "ERROR: cannot find applet version %s in the log" % ( type(dx_app_ver) )
            sys.exit(0)
        ana_step = self.enc_analysis_step_find(dx_app_name,dx_app_ver,dx_app_id,wf_run_notes["pipeline_name"])
        # step_run['analysis_step'] will be object coming from enc.
        if 'analysis_step' not in step_run:
            print "  + Need to update missing step_run['step'] new: '"+ana_step['name']+"'"
            step_run['analysis_step'] = "/analysis-steps/" + ana_step['name']
            update_step_run['analysis_step'] = "/analysis-steps/" + ana_step['name']
            patch_required = True
        else:
            step = step_run['analysis_step']
            if 'name' not in step or step['name'] != ana_step['name']:
                print "  + Need to update step_run['step'] new: '"+ana_step['name']+"'"
                #if 'name' in step:
                #    print "    enc: '"+step['name']+"'"
                #else:
                #    print json.dumps(step,indent=4,sort_keys=True)
                step_run['analysis_step'] = "/analysis-steps/" + ana_step['name']
                update_step_run['analysis_step'] = "/analysis-steps/" + ana_step['name']
                patch_required = True
            
        # Applet details are... detailed
        applet_details = {}
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
            
        # Applet details: layers of checking    
        if "dx_applet_details" not in step_run or len(step_run["dx_applet_details"]) == 0:
            print "  + Need to update missing step_run['dx_applet_details']"
            step_run["dx_applet_details"] = [ applet_details ]
            update_step_run["dx_applet_details"] =  [ applet_details ]
            patch_required = True
        else:
            assert len(step_run["dx_applet_details"]) == 1
            enc_applet_details = step_run["dx_applet_details"][0]
            for key in applet_details.keys():
                if key not in enc_applet_details:
                    print "  + Need to update missing step_run['dx_applet_details']['"+key+"'] new: '"+applet_details[key]+"'"
                    #print json.dumps(step_run["dx_applet_details"],indent=4,sort_keys=True)
                    step_run["dx_applet_details"] = [ applet_details ]
                    update_step_run["dx_applet_details"] =  [ applet_details ]
                    patch_required = True
                    break
                else:
                    if key != "parameters":
                        if applet_details[key] != enc_applet_details[key]:
                            print "  + Need to update step_run['dx_applet_details']['"+key+"']"
                            step_run["dx_applet_details"] = [ applet_details ]
                            update_step_run["dx_applet_details"] =  [ applet_details ]
                            patch_required = True
                            break
                    else:
                        param_changed = False
                        for param_key in applet_details[key].keys():
                            if param_key not in enc_applet_details[key] \
                            or applet_details[key][param_key] != enc_applet_details[key][param_key]:
                                param_changed = True
                                print "  + Need to update step_run['dx_applet_details']['"+key+"']['"+param_key+"']"
                                print json.dumps(step_run["dx_applet_details"],indent=4,sort_keys=True)
                                break
                        if param_changed:
                            step_run["dx_applet_details"] = [ applet_details ]
                            update_step_run["dx_applet_details"] =  [ applet_details ]
                            patch_required = True
                            break
            
        # Only updates notes if other things are being patched, and then always update!
        if patch_required:  
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
            update_step_run["notes"] = step_run["notes"]
            
            # Now patch this object
            self.obj_cache["exp"][step_alias] = step_run
            if test:
                print "  * Would patch step_run: '%s'" % step_alias
            else:
                try:
                    patched_step_run = dxencode.encoded_patch_obj(step_run['@id'], update_step_run, \
                                                                            self.server, self.authid, self.authpw)
                except:
                    print "Failed to patch step_run: '%s'" % step_alias
                    print json.dumps(update_step_run,indent=4,sort_keys=True)
                    sys.exit(1)
                print "  * Patched step_run: '%s'" % step_alias
                self.step_runs_patched += 1

        if step_run and verbose:
            sys.stderr.write("ENCODEd 'analysis_step_run':\n")
            sys.stderr.write(json.dumps(step_run,indent=4) + '\n')
            if "notes" in step_run:
                step_run_notes = json.loads(step_run.get("notes"))
                sys.stderr.write("ENCODEd 'analysis_step_run[notes]':\n")
                sys.stderr.write(json.dumps(step_run_notes,indent=4) + '\n')
        return step_run
        
    def file_metadata_recovery(self,fid,payload,test=True,verbose=False):
        '''Compares DX and encoded metadata and updates ENCODEd is necessary.'''
        recovered = False
        if fid not in self.found:
            print "* ERROR: Expecting to have exp_file for " + dxencode.file_path_from_fid(fid,projectToo=True)
            sys.exit(1)
            
        enc_file = self.found[fid]
        accession = enc_file['accession']
        update_payload = {}
        patch_required = False
        
        # Compare derived from and update if necessary
        derived_diffs = len(payload['derived_from'])
        if verbose:
            sys.stderr.write("> DX file derived_from:" + '\n')
            sys.stderr.write(json.dumps(payload['derived_from'],indent=4) + '\n')
        enc_derived = enc_file.get('derived_from')
        if enc_derived != None:
            if verbose:
                sys.stderr.write("> Enc file derived_from:" + '\n')
                sys.stderr.write(json.dumps(enc_derived,indent=4) + '\n')
            if len(enc_derived) == derived_diffs:
                for acc in payload['derived_from']:
                    for inp_file in enc_derived:
                        if acc == inp_file.get('accession'):
                            derived_diffs -= 1
        else:
            if verbose:
                sys.stderr.write("Enc file derived_from: Not Found\n")
        if derived_diffs > 0:
            print "  + Need to update 'derived_from'."
            update_payload['derived_from'] = payload['derived_from']
            patch_required = True
        
        # Compare wf_run/step_run and update if necessary
        step_run_diff = True
        dx_step_run = payload.get("step_run")
        if dx_step_run != None:
            dx_step_run = dx_step_run.split('/')[-1] # Just the alias please
            if verbose:
                sys.stderr.write("> DX file step_run:" + dx_step_run  + '\n')
        else:
            print "* ERROR: payload is missing 'step_run'"
        enc_step_run = enc_file.get("step_run")
        if verbose:
            if enc_step_run != None:
                sys.stderr.write("> Enc file step_run:" + enc_step_run +'\n')
                #sys.stderr.write(json.dumps(enc_step_run,indent=4) + '\n')
            else:
                sys.stderr.write("> Enc file step_run: Not Found\n")

        if dx_step_run == None and enc_step_run == None:   
            step_run_diff = False
        elif dx_step_run != None and enc_step_run != None:
            step_run = self.enc_lookup_json(enc_step_run,must_find=True)
            if verbose:
                sys.stderr.write("> Actual step_run:" + '\n')
                sys.stderr.write(json.dumps(step_run,indent=4) + '\n')
            aliases = step_run.get('aliases')
            if aliases != None:
                if verbose:
                    sys.stderr.write("> aliass:" + '\n')
                    sys.stderr.write(json.dumps(aliases,indent=4) + '\n')
                if dx_step_run in aliases:
                    step_run_diff = False
        if step_run_diff:
            print "  + Need to update 'step_run'."
            update_payload["step_run"] = payload["step_run"]
            patch_required = True
            
        # What about lab?
        if 'lab' not in enc_file:
            print "  + Need to update missing 'lab' new: '"+payload['lab']+"'"
            update_payload['lab'] = payload['lab']
            patch_required = True
        else:
            lab = enc_file.get('lab')
            if '@id' not in lab or lab['@id'] != payload['lab']:
                print "  + Need to update 'lab' new: '"+payload['lab']+"'"
                print json.dumps(lab,indent=4,sort_keys=True)
                update_payload['lab'] = payload['lab']
                patch_required = True

        # award?
        if 'award' not in enc_file:
            print "  + Need to update missing 'award' new: '"+payload['award']+"'"
            update_payload['award'] = payload['award']
            patch_required = True
        elif enc_file['award'] != payload['award']:
            print "  + Need to update 'award' new: '"+payload['award']+"'  enc: '"+enc_file['award']+"'"
            update_payload['award'] = payload['award']
            patch_required = True

        if patch_required:
            update_payload['notes'] = payload['notes'] # Only update notes if other things update, then always update notes
            if test:
                print "  * Would patched file: '%s'" % dxencode.file_path_from_fid(fid)
                #print json.dumps(update_payload,indent=4,sort_keys=True)
                recovered = False  # would recover
            else:
                try:
                    ret = dxencode.encoded_patch_obj(accession, update_payload, self.server, self.authid, self.authpw)
                except:
                    print "Failed to patch file: '%s'" % dxencode.file_path_from_fid(fid)
                    print json.dumps(update_payload,indent=4,sort_keys=True)
                    sys.exit(1)
                print "  * Patched file: '%s'" % dxencode.file_path_from_fid(fid)
                recovered = True
                
        return recovered

    def run(self):
        '''Runs recovery from start to finish using command line arguments.'''
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
        print "== Running in project [%s] and will attempt recovery to the [%s] server ==" % \
                                                        (self.proj_name,self.server_key)

        exp_count = 0
        halted = 0
        total_recovered = 0
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
            print "- Found %d files that are available in DX." % len(files_expected)
            if len(files_expected) == 0:
                continue

            # 5) For each file that should be posted, determine if the file needs to be posted.
            files_posted = self.find_posted_files(files_expected, verbose=args.verbose) #True)
            print "- Found %d files that have been posted" % len(files_posted)
            if len(files_posted) == 0:
                continue

            # 6) For each file that needs to be posted:
            exp_count += 1
            file_count = 0
            recovery_count = 0
            for (out_type,rep_tech,fid) in files_posted:
                sys.stdout.flush() # Slow running job should flush to piped log
                # a) discover all necessary dx information needed for post.
                # b) gather any other information necessary from dx and encoded.
                accession = self.found[fid]['accession']
                if args.start_at != None:
                    if accession != args.start_at:
                        continue
                    else:
                        print "- Starting at %s" % (accession)
                        args.start_at = None
                    
                print "- Handle file %s %s" % (accession,dxencode.file_path_from_fid(fid))
                payload = self.make_payload_obj(out_type,rep_tech,fid, verbose=args.verbose)

                file_count += 1
                # c) Update encoded database only if necessary.
                if self.file_metadata_recovery(fid,payload,args.test,verbose=args.verbose):
                    if not args.test:
                        recovery_count += 1

                if args.files != 0 and file_count >= args.files:  # Short circuit for test
                    print "- Just trying %d file(s) by request" % file_count
                    break

            print "- For %s Processed %d file(s), recovered %s" % \
                                                        (self.exp_id, file_count, recovery_count)
            total_recovered += recovery_count

        print "Processed %d experiment(s), halted %d, recovered %d file(s)" % \
                                                            (exp_count, halted, total_recovered)
        if halted == exp_count:
            sys.exit(1)
        print "(finished)"


if __name__ == '__main__':
    '''Run from the command line.'''
    recover = Recovery()
    recover.run()

