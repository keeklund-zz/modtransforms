from os.path import exists

def update_header(header, chrom_sizes_file):
    """Return updates header based on chromosome sizes file.

    Function takes header data and returns a modified data structure with 
    values taken from chromosome sizes file.

    """
    assert isinstance(header, dict), "header structure must be type dict"
    assert exists(chrom_sizes_file), "chromosome sizes file must exist"
    with open(chrom_sizes_file, "r") as cs_file:
        chrom_sizes = dict(line.split() for line in cs_file)
    for chrom in header['SQ']:
        chrom['LN'] = int(chrom_sizes[chrom['SN']])
    return header
