import gzip

def gen_mod(mod_file):
    """Return mod file generator"""
    with gzip.open(mod_file, 'r') as mod:
        for line in mod:
            if not line.startswith('#'):
                yield line.rstrip()

