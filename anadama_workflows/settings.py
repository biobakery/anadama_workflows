class workflows:
    product_directory = "mibc_products"
    class metaphlan2(object):
        bowtie2db = "/seq/ibdmdb/centos6/var/lib/bowtie2db/mpa.200.ffn"
        mpa_pkl   = "/seq/ibdmdb/centos6/var/lib/bowtie2db/mpa.200.pkl"
    class sixteen(object):
        otu_taxonomy = "/seq/ibdmdb/centos6/var/lib/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt"
        otu_refseq   = "/seq/ibdmdb/centos6/var/lib/gg_13_5_otus/rep_set/97_otus.fasta"
    class alignment:
        kegg_bowtie2_db = "/vagrant/data/ccfa/KEGG_reduced_bowtie2/KEGG_1"

