#!/usr/bin/env python

import os
import sys
import dxencode
import argparse

class Checker(object):
    ''' Abstract class, subclass for each assay/pipeline '''
    def __init__(self):
        self.args = self.get_args()
        # expect graph of form (child, parent)
        self.expected_graph = {
            'unpaired': (()),
            'paired': (())
        }
        self.assay_term_name = ''
        if args.test:
            key = 'test'
        else:
            key = 'www'
        (self.authid, self.authpw, self.server) = dxencode.processkey(key)

    def run(self):
        self.experiments = self.do_query(self.args)
        for exp in self.experiments:
            self.reconcile(exp)

    def do_query(self, args):
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
                        required=False)

        ap.add_argument('-a', '--all',
                        help='All experiments for assay type (as specified in self.do_query())',
                        nargs='+',
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

    def reconcile(self, exp):

        reads = {}
        derived = {}
        paired = 'unpaired'
        if dxencode.is_paired_ended(exp):
            paired = 'unpaired'

        for f in exp['files']:  ## this could throw a key error but shouldn't
            if f['file_format'] == 'fastq' and f['output_type'] == 'reads':
                rstr = 'rep_'+f['biological_replicate_number']+'_'+f['technical_replicate_number']
                rr = reads.get(rstr,{}).update({f['accession']: f})
                reads[rstr] = rr
            elif f.get('notes', ''):
                try:
                    notes = json.loads(f['notes'])
                    # this should only match pipeline generated files
                    try:
                        gen = f['assembly']
                    except:
                        print("File has no assembly: %s %s %s %s %s" % (f['accession'], f['file_format'], f['output_type'], f['submitted_by']['email'], notes))
                        continue
                    try:
                        ann = f['genome_annotation']
                    except:
                        print("File has no annotation: %s %s %s %s %s" % (f['accession'], f['file_format'], f['output_type'], f['submitted_by']['email'], notes))
                        continue

                    rfg = derived.get(gen,{})
                    rfa = rfg.get(ann, {}).update({f['accession']: f})

                    derived[gen][ann] = rfa
                    # do I need to hash by output_type here?
                except:
                    print("Other file: %s %s %s %s" % (f['accession'], f['file_format'], f['output_type'], f['submitted_by']['email']))
                    continue

        for assembly in derived.keys():
            for annotation in genome.keys():
                for f in derived[assembly][annotation].values():
                    dfs = f['derived_from']
                    if not dfs:
                        print("File not in graph: %s %s %s %s %s" % (f['accession'], f['file_format'], f['output_type'], f['submitted_by']['email'], notes))
                        print("Possible parents: %s" %(self.guess_parents(f, derived[assembly][annotation], reads)))
                    else:
                        current_graph = self.expected_graph[paired]
                        notes = json.loads(f['notes'])
                        f_rstr = 'rep_'+f['biological_replicate_number']+'_'+f['technical_replicate_number']
                        aligns = [ sv['software'] for sv in notes['software_versions'] if sv['software'] in self.aligners ]

                        if f['output_type'] == 'alignments':
                            if len(aligns) == 1:
                                aligner = aligner[0]
                            else:
                                print("File %s has too many aligners %s" % (f['accession'], notes['software_versions']))
                                continue
                        for df in dfs:
                            try:
                                assert(df['output_type'] == current_graph[f['output_type']])
                                del current_graph[f['output_type']]
                                # only one entry per genome/assembly
                            except AssertionError:
                                print("File %s (%s) has mismatched output types with derived from %s (%s)" % (f['accession'], f['output_type'], df['accession'], df['output_type']))

                            ## Note this will not be true of all files in assays
                            try:
                                df_rstr = 'rep_'+df['biological_replicate_number']+'_'+df['technical_replicate_number']
                                assert(df_rstr == f_rstr)
                            except AssertionError:
                                print("File %s (%s) is derived from file from different replicate %s (%s)" %(f['accesssion'], f_rstr, df['accession'], df_rstr))

                            if df['file_format'] == 'bigWig':
                                try:
                                    assert(df['submitted_file_name'].find(aligner) >= 0) # God have mercy on my soul
                                except AssertionError:
                                    print("Alignement File %s (%s) doesn't match it's bigWig %s" %(f['accesssion'],aligner,df['submitted_file_name']))



    def do_query(self, args):
        if args.experiments:
            self.experiments = []
            for acc in args.experiments:
                r = dxencode.encoded_get(self.server+acc, AUTHID=self.authid,AUTHPW=self.authpw)
                try:
                    r.raise_for_status()
                    e = r.json()
                    if e.get('replicates',[]) and [ f for f in e.get('files',[]) if f['file_format'] == 'fastq' ]:
                        self.experiments.append(e)
                except:
                    print("Could not find %s in encodeD" % (acc))
                    continue
        elif args.all:
            q = self.server+'search/?type=experiment&assay_term_id=%s&award.rfa=ENCODE3&limit=all&files.file_format=fastq' % self.ASSAY_TERM_ID
            r = dxencode.encoded_get(self.server+q,  AUTHID=self.authid,AUTHPW=self.authpw)
            try:
                r.raise_for_status()
                self.experiments = [ e for e in r.json()['@graph'] if e['replicates'][0]['library'].get('size_range', "") != '>200' ]
            except Exception, ex:
                print("Some error: %s trying to find experiments" % (ex))
                raise

    def __init__(self):
        self.expected_graph = {
        'unpaired':
            [
                {'alignments': 'reads'},
                {'transcriptome alignments': 'reads'},
                {'genome quantifications': 'transcriptome alignments'},
                {'transcript quantifications': 'transcriptome alignments'},
                {'unique signal': 'alignments'},
                {'multi-read signal': 'alignments'}
            ],
        'paired':
            [
                {'alignments': 'reads'},
                {'transcriptome alignments': 'reads'},
                {'genome quantifications': 'transcriptome alignments'},
                {'transcript quantifications': 'transcriptome alignments'},
                {'unique plus signal': 'alignments'},
                {'multi-read plus signal': 'alignments'},
                {'unique minus signal': 'alignments'},
                {'multi-read minus signal': 'alignments'}
           ],
        }
        self.aligners = ['TopHat', 'STAR']

if __name__ == '__main__':
    '''Run from the command line.'''
    check = lRNAChecker()
    check.run()
