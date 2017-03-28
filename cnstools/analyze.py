import os, json
import file_handlers as fhs

def analyze(out, analyses, align=None,score=None,identify=None):
    analyzer = _Analyzer(out,align,score,identify)
    for analysis in analyses:
        getattr(analyzer,analysis)()
    print json.dumps(analyzer.results,sort_keys=True,indent=4)

class _Analyzer(object):
    """docstring for _Analyzer"""
    def __init__(self, out, align=None,score=None,identify=None):
        super(_Analyzer, self).__init__()
        self.original_wd = os.getcwd()
        self.out = out
        self.align = align
        self.score = score
        self.identify = identify
        self.results = {}
    def align_coverage(self):
        if not self.align:
            raise ValueError("Did not specify alignment results file.")
        coverage = {}
        for genome in self.align["aligned_query_genomes"]:
            r_srcs = set()
            q_srcs = set()
            ref_size = 0
            que_size = 0
            ref_bases, que_bases = [0]*2
            for path in self.align["aligned_query_genomes"][genome]:
                alignment = fhs.maf.Handler(path)
                for entry in alignment.entries(include_comments=False):
                    if entry.sequences[1].src not in q_srcs:
                        que_size += entry.sequences[1].srcSize
                        q_srcs.add(entry.sequences[1].src)
                    if entry.sequences[0].src not in r_srcs:
                        ref_size += entry.sequences[0].srcSize
                        r_srcs.add(entry.sequences[0].src)
                    ref_bases += entry.sequences[0].size
                    que_bases += entry.sequences[1].size
            coverage[genome] = {"reference":ref_bases/float(ref_size),"query":que_bases/float(que_size)}
        self.results["align_coverage"] = coverage
            





        

def cl_analyze(out, analyses, align_results=None,score_results=None,identify_results=None):
    print align_results, out, analyses
    align, score, identify = [None]*3
    if align_results!=None:
        with open(align_results) as align_file:
            align = json.loads(align_file.read())
            align["dirname"] = os.path.dirname(os.path.abspath(align_results))
    if score_results!=None:
        with open(score_results) as score_file:
            score = json.loads(score_file.read())
            score["dirname"] = os.path.dirname(os.path.abspath(score_results))
    if identify_results!=None:
        with open(identify_results) as identify_file:
            identify = json.loads(identify_file.read())
            identify["dirname"] = os.path.dirname(os.path.abspath(identify_results))
    analyze(out, analyses, align, score, identify)

_cl_entry = cl_analyze #function that should be run on command line entry to this subcommand
def _parser(parser_add_func,name):
    p = parser_add_func(name,description="Aligns a genome to a reference")
    fg = p.add_mutually_exclusive_group(required=True)
    fg.add_argument("-a","--align_results", help="`align` results JSON file")
    fg.add_argument("-s","--score_results", help="`align` results JSON file")
    fg.add_argument("-i","--identify_results", help="`align` results JSON file")
    p.add_argument("-o", "--out", required=True, help="output folder")
    p.add_argument("analyses", help="Which analyses to run",nargs="+")
    return p

    # #for iteration over genomes regardless of reference
    # all_genomes = config["genomes"]

    # #maf_to_bed
    # info = "Convert per-genome CNS regions to .bed:"
    # header_print(info)
    # chrom['genome_cns_beds_folder'] = os.path.join(chrom_out,"genome_cns_beds")
    # try:
    #     os.makedirs(chrom['genome_cns_beds_folder'])
    # except OSError:
    #     if not os.path.isdir(chrom['genome_cns_beds_folder']):
    #         raise
    # chrom['genome_cns_beds'] = {}
    # cns_maf = Maf(file_name=chrom['cns_maf'])
    # for genome in all_genomes:
    #     chrom['genome_cns_beds'][genome] = os.path.join(chrom['genome_cns_beds_folder'],genome+"_cns_"+chrom_name+".bed")
    #     bed = cns_maf.to_bed(genome_name=genome,index_tag="cns_maf_index")
    #     bed.save_file(chrom['genome_cns_beds'][genome])
    # del cns_maf


    # #$bedtools closest
    # info = "Find closest gene for each CNS region:"
    # header_print(info)
    # chrom['gene_proximity_beds_folder'] = os.path.join(chrom_out,"gene_proximity_beds")
    # try:
    #     os.makedirs(chrom['gene_proximity_beds_folder'])
    # except OSError:
    #     if not os.path.isdir(chrom['gene_proximity_beds_folder']):
    #         raise
    # chrom['gene_proximity_beds'] = {}
    # for genome in all_genomes:
    #     chrom['gene_proximity_beds'][genome] = os.path.join(chrom['gene_proximity_beds_folder'],genome+"_proxim.bed")
    #     cmd = "bedtools closest -D a -a %s -b %s > %s" % \
    #         (chrom['genome_cns_beds'][genome],
    #          config['genome_gene_beds'][genome],
    #          chrom['gene_proximity_beds'][genome])
    #     process = subprocess.Popen(cmd, shell=True)
    #     process.wait()

    # #maf_and_proxim_bed_to_cns
    # info = "Process proximity and maf files into .cns file:"
    # header_print(info)
    # chrom['results'] = os.path.join(chrom_out,"identified_CNSs.cns")
    # cns_proxim_beds = {genome:Bed13(chrom['gene_proximity_beds'][genome]) for genome in all_genomes}
    # Maf(file_name=chrom['cns_maf'])\
    #     .cns_from_proxim_beds(cns_proxim_beds,"cns_maf_index",prefix=chrom_name+".")\
    #     .save_file(chrom['results'])   