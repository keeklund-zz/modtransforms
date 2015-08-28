from modtransforms.formats.process import bam, gtf, npf, smrna_12_bed, smrna_bed
from modtransforms.formats.process import smrna_gff3, smrna_lib_fa
from modtransforms.formats.process import smrna_table_txt, smrna_txt

file_types = {
    'bam': bam,
    'gtf': gtf,
    'npf': npf,
    'smrna_12_bed': smrna_12_bed,
    'smrna_bed': smrna_bed,
    'smrna_gff3': smrna_gff3,
    'smrna_lib_fa': smrna_lib_fa,
    'smrna_table_txt': smrna_table_txt,
    'smrna_txt': smrna_txt,
}

