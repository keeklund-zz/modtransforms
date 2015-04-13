import gzip

def gen_mod(mod_file):
    """Return mod file generator"""
    with gzip.open(mod_file, 'r') as mod:
        for line in mod:
            if not line.startswith('#'):
                yield line.rstrip()

def build_transform(mod_file):
    """Return nested dictionary: chr: pos1: pos2"""
    mod = gen_mod(mod_file)
    transform = {}

    chrom = ''
    while mod:
        try:
            data = mod.next().split()
        except:
            break
        if data[1] != chrom:
            delta = 0 
            chrom = data[1]

        if data[0] == 's':
            pass
        elif data[0] == 'd':
            delta = delta + -1
        elif data[0] == 'i':
            delta = delta + len(data[3])
        else:
            exit('not controlling every option')

        try:
            transform[chrom].append((int(data[2]), delta))
        except:
            transform[chrom] = [(int(data[2]), delta),]
    return transform

