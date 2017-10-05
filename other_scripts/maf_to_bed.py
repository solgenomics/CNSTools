import cnstools.file_handlers.maf as maf
import sys

maf.Handler(sys.argv[1]).to_bed(sys.argv[2],
                                 genome_name="Alyra",
                                 index_tag="cns.maf_index",
                                 tracker_name="Converting")