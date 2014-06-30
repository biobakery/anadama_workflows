class workflows(object):
    product_directory = "mibc_products"

    class metaphlan2(object):
        bowtie2db = "/vagrant/metaphlandb/mpa.200.ffn"
        mpa_pkl   = "/vagrant/metaphlandb/mpa.200.pkl"
    class sixteen(object):
        otu_taxonomy = "/vagrant/data/ccfa/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt"
        otu_refseq   = "/vagrant/data/ccfa/gg_13_8_otus/rep_set/99_otus.fasta"
    class alignment(object):
        kegg_bowtie2_db = "/vagrant/data/ccfa/KEGG_reduced_bowtie2/KEGG_1"


