#!/usr/bin/env python
import splashdown
import dxencode
import sys
import dxpy.exceptions

class Patchdown(splashdown.Splashdown):

    post_templates = {
        # For looking up previous result files, use wild-cards
        # used with old lrnaLaunch script
        "tophat_bam":   {
            "file_format": "bam",
            "output_type": "alignments",
            "derived_from": ["reads1", "reads2"]
        },
        "tophat_minus_all_bw":  {
            "file_format": "bigWig",
            "output_type": "multi-read minus signal",
            "derived_from": ["tophat_bam"]
        },
        "tophat_minus_uniq_bw":{
            "file_format": "bigWig",
            "output_type": "unique minus signal",
            "derived_from": ["tophat_bam"]
        },
        "tophat_plus_all_bw":   {
            "file_format": "bigWig",
            "output_type": "multi-read plus signal",
            "derived_from": ["tophat_bam"]
        },
        "tophat_plus_uniq_bw":  {
            "file_format": "bigWig",
            "output_type": "unique minus signal",
            "derived_from": ["tophat_bam"]
        },
        "tophat_all_bw":        {
            "file_format": "bigWig",
            "output_type": "multi-read signal",
            "derived_from": ["tophat_bam"]
        },
        "tophat_uniq_bw":       {
            "file_format": "bigWig",
            "output_type": "unique signal",
            "derived_from": ["tophat_bam"]
        },
        "star_genome_bam":      {
            "file_format": "bam",
            "output_type": "alignments",
            "derived_from": ["reads1", "reads2"]
        },
        "star_minus_all_bw":    {
            "file_format": "bigWig",
            "output_type": "multi-read minus signal",
            "derived_from": ["star_genome_bam"]
        },
        "star_minus_uniq_bw":   {
            "file_format": "bigWig",
            "output_type": "unique minus signal",
            "derived_from": ["star_genome_bam"]
        },
        "star_plus_all_bw":     {
            "file_format": "bigWig",
            "output_type": "multi-read plus signal",
            "derived_from": ["star_genome_bam"]
        },
        "star_plus_uniq_bw":    {
            "file_format": "bigWig",
            "output_type": "unique plus signal",
            "derived_from": ["star_genome_bam"]
        },
        "star_all_bw":          {
            "file_format": "bigWig",
            "output_type": "multi-read signal",
            "derived_from": ["star_genome_bam"]
        },
        "star_uniq_bw":         {
            "file_format": "bigWig",
            "output_type": "unique signal",
            "derived_from": ["star_genome_bam"]
        },
        "rsem_gene_results":    {
            "file_format": "tsv",
            "output_type": "genome quantifications",
            "derived_from": ["star_anno_bam"]
            # note should be derived from star_anno_bam
        },
        "star_anno_bam":        {
            "file_format": "bam",
            "output_type": "transcriptome alignments",
            "derived_from": ["reads1", "reads2"]
        },
        "rsem_iso_results":     {
            "file_format": "tsv",
            "output_type": "transcript quantifications",
            "derived_from": ["star_anno_bam"]
        },
        "reads1":     {
            "file_format": "fastq",
            "output_type": "reads1",
            "derived_from": []
        },
        "reads2":     {
            "file_format": "fastq",
            "output_type": "reads2",
            "derived_from": []
        }
    }

    def __init__(self):
        super(Patchdown, self).__init__()
        self.derived_map = { x['output_type']:
            [ self.post_templates[y]['output_type'] for y in x['derived_from'] ] for x in self.post_templates.values() }


    def find_derived_from(self,fid,job,verbose=False):
        ''' wrap in exception because sometimes dx-files go missing '''
        try:
            derived_from = super(Patchdown, self).find_derived_from(fid,job,verbose)

        except dxpy.exceptions.ResourceNotFound, e:
            print "WARN: derived_from failed %s" % e
            derived_from = []

        if not derived_from:
            #import pdb;pdb.set_trace()
            # try to guess
            pass
        return derived_from

    def run(self):
        '''Override super.run()'''
        args = self.get_args()
        self.test = args.test
        self.server_key = args.server
        if self.server_key != "test":
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
        for exp_id in args.experiments:
            sys.stdout.flush() # Slow running job should flush to piped log
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
            self.exp_folder = dxencode.find_exp_folder(self.project,exp_id,args.results_folder,warn=True)
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
            files_expected = self.find_expected_files(self.exp_folder,self.replicates, verbose=args.verbose)
            print "- Found %d files that are available to post." % len(files_expected)
            if len(files_expected) == 0:
                continue

            # 5) For each file that should be posted, determine if the file needs to be posted.
            files_to_post = { x[2]: x for x in self.find_needed_files(files_expected, verbose=args.verbose) }
            # index on dx file id
            print "- Found %d files that need to be posted" % len(files_to_post.keys())

            # 6) For each file that needs to be posted:
            exp_count += 1
            file_count = 0
            post_count = 0
            for (out_type,rep_tech,fid) in files_expected:
                sys.stdout.flush() # Slow running job should flush to piped log
                # a) discover all necessary dx information needed for post.
                # b) gather any other information necessary from dx and encoded.
                print "  Handle file %s" % dxencode.file_path_from_fid(fid)
                job = dxencode.job_from_fid(fid)

                try:
                    derived_from = self.find_derived_from(fid,job, args.verbose)

                except dxpy.exceptions.ResourceNotFound, e:
                    print "WARN: derived_from failed %s" % e
                    derived_from = []
                if not files_to_post.get(fid,()):
                    f_obj = self.found.get(fid,None)
                    if f_obj:
                        current_derived_from = f_obj['derived_from']
                        if derived_from and not current_derived_from:
                            print "Need to patch derived_from for %s/%s to %s (currently: %s)" % (f_obj['accession'], fid, derived_from, current_derived_from)
                        else:
                            print "Derived from for %s good" % f_obj['accession']
                    else:
                        print "File %s (%s) from %s/%s not found @ DNANexus" % (fid,out_type,exp_id,rep_tech)

                #POSTING
                else:
                    payload = self.make_payload_obj(out_type,rep_tech,fid, verbose=args.verbose)
                    if args.force_annotation:
                        print "WARN: forcing genome_annotation to be %s" % args.force_annotation
                        payload['genome_annotation'] = args.force_annotation

                    file_count += 1
                    # c) Post file and update encoded database.
                    accession = self.file_post(fid,payload,args.test)
                    if accession == None:
                        print "* HALTING %s - post failure could compromise 'derived_from'" % \
                                                                                        (self.exp_id)
                        halted += 1
                        break

                    # d) Update dnanexus file with file accession tag.
                    if not args.test:
                        post_count += 1
                    self.file_mark_accession(fid,accession,args.test)

                print "- For %s Processed %d file(s), posted %s" % \
                                                            (self.exp_id, file_count, post_count)
                total_posted += post_count
        print "Processed %d experiment(s), halted %d, posted %d file(s)" % \
                                                            (exp_count, halted, total_posted)
        if halted == exp_count:
            sys.exit(1)
        print "(finished)"


if __name__ == '__main__':
    '''Run from the command line.'''
    patch = Patchdown()
    patch.run()


