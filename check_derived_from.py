#!/usr/bin/env python

import os
import sys
import json
import dxencode
import argparse

class Checker(object):
    ''' Abstract class, subclass for each assay/pipeline '''
    def __init__(self):
        self.args = self.get_args()
        # expect graph of form (child, parent)
        self.expected_graph = {
            'unpaired': {},
            'paired': {}
        }
        self.assay_term_name = ''
        if self.args.test:
            key = 'test'
        else:
            key = 'www'
        (self.authid, self.authpw, self.server) = dxencode.processkey(key)
        self.experiments = []

    def run(self):
        self.do_query()
        if not self.experiments:
            print "No experiments found."
            sys.exit(1)
        for exp in self.experiments:
            self.reconcile(exp)

    def do_query(self):
        ''' returns a list of encodeD experiment objects given cmd line args'''
        pass

    def reconcile(self, exp):
        ''' compares existing derived from graph(s) with ideal one specified in subclass '''
        pass

    def get_args(self):
        '''Parse the input arguments.'''
        ### PIPELINE SPECIFIC
        ap = argparse.ArgumentParser(description="For a given experiment or set of experiments loops over the pipeline created files and sees if they are consistant")
        ### PIPELINE SPECIFIC

        ap.add_argument('-e', '--experiments',
                        help='One or more ENCODED experiment accessions',
                        nargs='+',
                        default=[],
                        required=False)

        ap.add_argument('-a', '--all',
                        help='All experiments for assay type (as specified in self.do_query())',
                        action='store_true',
                        required=False)

        ap.add_argument('--verbose',
                        help='More debugging output.',
                        action='store_true',
                        required=False)

        ap.add_argument('--test',
                        help='More debugging output.',
                        action='store_true',
                        required=False)

        return ap.parse_args()

class lRNAChecker(Checker):
    ASSAY_TYPE = 'RNA-seq'
    ASSAY_TERM_ID = "OBI:0001271"

    def guess_parents(self, f, derived_files, reads, current_graph):
        parent_type = current_graph[f['output_type']]
        possibles = [ (d['accession'], d['file_format'], d['output_type'], d['submitted_file_name']) for d in derived_files.values()
            if d['output_type'] == parent_type ]
        return possibles

    def reconcile(self, exp):

        reads = {}
        derived = {}

        for facc in exp.get('original_files',[]): ## could have no files
            res = dxencode.encoded_get(self.server+facc, AUTHID=self.authid, AUTHPW=self.authpw)
            try:
                res.raise_for_status()
                f = res.json()
            except:
                print("File: %s not found" % (facc))
                continue

            if f['file_format'] == 'fastq' and f['output_type'] == 'reads':
                try:
                    rstr = "rep_%s_%s" % (f['replicate']['biological_replicate_number'], f['replicate']['technical_replicate_number'])
                except:
                    print("missing replicate numbers for fastq: %s" % (json.dumps(f,indent=4)))
                    continue
                rr = reads.get(rstr,{})
                rr.update({f['accession']: f})
                reads[rstr] = rr
            else:
                try:
                    notes = json.loads(f['notes'])
                    # this should only match pipeline generated files
                    try:
                        gen = f['assembly']
                    except:
                        print("File has no assembly: %s %s %s %s %s" % (f['accession'], f['file_format'], f['output_type'], f['submitted_by']['title'], notes))
                        continue
                    try:
                        ann = f['genome_annotation']
                    except:
                        print("File has no annotation: %s %s %s %s %s" % (f['accession'], f['file_format'], f['output_type'], f['submitted_by']['title'], notes))
                        continue

                    rfg = derived.get(gen,{})
                    rfa = rfg.get(ann, {})
                    rfa.update({f['accession']: f})

                    if not derived.has_key(gen):
                        derived[gen] = {}

                    if derived[gen].has_key(ann):
                        derived[gen][ann] = rfa
                    else:
                        derived[gen].update({ ann: rfa})
                    # do I need to hash by output_type here?
                except (ValueError, KeyError), e:
                    # should also trap JSON loads exception.
                    print("WARN: Other file: %s %s %s %s" % (f['accession'], f['file_format'], f['output_type'], f['submitted_by']))
                    #print("Error: %s" % (e))
                    continue

        if not reads.keys():
            print("WARN: No read files found in %s" % (exp['accession']))
            return

        for rep in reads.keys():
            print("%s: Found %s fastqs in replicate %s" % (exp['accession'], len(reads[rep]), rep))


        if not derived.keys():
            print("WARN: No derived files found in %s" % (exp['accession']))
            return

        for gen in derived.keys():
            for ann in derived[gen].keys():
                print("%s: Found %s derived files for %s/%s" %(exp['accession'], len(derived[gen][ann]), gen,ann))

        paired = 'unpaired'
        if dxencode.is_paired_ended(exp):
            paired = 'paired'

        for assembly in derived.keys():
            for annotation in derived[assembly].keys():
                current_graph = self.expected_graph[paired]

                for acc, f in derived[assembly][annotation].items():
                    dfs = f['derived_from']
                    if not dfs:
                        print("File not in graph: %s %s %s %s %s" % (f['accession'], f['file_format'], f['output_type'], f['submitted_by']['title'], notes))
                        print("Possible parents: %s" %(self.guess_parents(f, derived[assembly][annotation], reads, current_graph)))
                    else:
                        notes = json.loads(f['notes'])
                        try:
                            f_rstr = "rep_%s_%s" % (f['replicate']['biological_replicate_number'], f['replicate']['technical_replicate_number'])
                        except:
                            print("missing replicate numbers for %s" % (json.dumps(f,indent=4)))
                            continue
                        aligns = [ sv['software'] for sv in notes['software_versions'] if sv['software'] in self.aligners ]

                        if f['output_type'] == 'alignments':
                            if len(aligns) == 1:
                                aligner = aligns[0]
                            else:
                                print("File %s has too many aligners %s" % (f['accession'], notes['software_versions']))
                                continue
                        for df in dfs:
                            try:
                                assert(df['output_type'] == current_graph[f['output_type']])
                                #del current_graph[f['output_type']]
                                # only one entry per genome/assembly
                            except AssertionError:
                                print("File %s (%s) has mismatched output types with derived from %s (%s)" % (f['accession'], f['output_type'], df['accession'], df['output_type']))

                            ## Note this will not be true of all files in assays
                            '''
                            ## Skipping this because derived from file doesn't have replicate embedded
                            try:
                                df_rstr = "rep_%s_%s" % (df['replicate']['biological_replicate_number'], df['replicate']['technical_replicate_number'])
                                assert(df_rstr == f_rstr)
                            except AssertionError:
                                print("File %s (%s) is derived from file from different replicate %s (%s)" %(f['accesssion'], f_rstr, df['accession'], df_rstr))
                            '''
                            if df['file_format'] == 'bigWig':
                                try:
                                    assert(df['submitted_file_name'].find(aligner) >= 0) # God have mercy on my soul
                                except AssertionError:
                                    print("Alignment File %s (%s) doesn't match it's bigWig %s" %(f['accesssion'],aligner,df['submitted_file_name']))



    def do_query(self):
        if self.args.experiments:
            self.experiments = []
            for acc in self.args.experiments:
                res = dxencode.encoded_get(self.server+acc, AUTHID=self.authid,AUTHPW=self.authpw)
                try:
                    #import pdb;pdb.set_trace()
                    res.raise_for_status()
                    e = res.json()
                    if e.get('replicates',[]) and [ f for f in e.get('files',[]) if f['file_format'] == 'fastq' ]:
                        self.experiments.append(e)
                except:
                    print("Could not find %s in encodeD" % (acc))
                    continue
        elif self.args.all:
            q = 'search/?type=experiment&assay_term_id=%s&award.rfa=ENCODE3&limit=all&files.file_format=fastq&frame=embedded' % self.ASSAY_TERM_ID
            res = dxencode.encoded_get(self.server+q,  AUTHID=self.authid,AUTHPW=self.authpw)
            try:
                res.raise_for_status()
                self.experiments = [ e for e in res.json()['@graph'] if e['replicates'][0]['library'].get('size_range', "") != '>200' ]
            except Exception, ex:
                print("Some error: %s trying to find experiments" % (ex))
                raise
        else:
            print "Specify -e/--experiments for a list of experiments or --all for all long RNA seq expts"
            sys.exit(1)

        print("%s Experiments" %(len(self.experiments)))

    def __init__(self):
        super(lRNAChecker, self).__init__()

        self.expected_graph = {
        'unpaired':
            {
                'alignments': 'reads',
                'transcriptome alignments': 'reads',
                'genome quantifications': 'transcriptome alignments',
                'transcript quantifications': 'transcriptome alignments',
                'unique signal': 'alignments',
                'multi-read signal': 'alignments'
            },
        'paired':
                {
                    'alignments': 'reads',
                    'transcriptome alignments': 'reads',
                    'genome quantifications': 'transcriptome alignments',
                    'transcript quantifications': 'transcriptome alignments',
                    'unique plus signal': 'alignments',
                    'multi-read plus signal': 'alignments',
                    'unique minus signal': 'alignments',
                    'multi-read minus signal': 'alignments'
                }
        }
        self.aligners = ['TopHat', 'STAR']


if __name__ == '__main__':
    '''Run from the command line.'''
    check = lRNAChecker()
    check.run()
