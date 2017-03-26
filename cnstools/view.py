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