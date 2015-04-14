from pysam import AlignmentFile
from modtransforms.mod.transform import build_transform, find_delta

def gtf(args, logger):
    out = open(args.output, 'w')
    chrom_mods = build_transform(args.mod, logger, args.reverse)
    keys = 'contig feature source start end score strand frame attributes'
    curr_chrom = ""
    with open(args.input, 'r') as input_:
        for line in input_:
            gene = line.rstrip()
            data = {k:g for k,g in zip(keys.split(), gene.split('\t'))}
            if data.get('contig') != curr_chrom:
                curr_chrom = data.get('contig')
                try:
                    positions, deltas = zip(*chrom_mods.get(curr_chrom))
                    logger.info("CONTIG: '%s'" % curr_chrom)
                except TypeError:
                    logger.warn(
                        "CONTIG: '%s' is not in MOD File. Skipping." % \
                        curr_chrom)
                    positions, deltas = [], 0
                    continue
            try:
                start_delta = find_delta(positions, 
                                         deltas, 
                                         int(data.get('start')))
                end_delta = find_delta(positions, 
                                       deltas,
                                       int(data.get('end')))
                data['start'] = int(data.get('start')) + start_delta
                data['end'] = int(data.get('end')) + end_delta
                out.write('%s\n' % \
                          '\t'.join([str(data.get(k)) for k in keys.split()]))
            except IndexError:
                pass
            
def smrna_gff3(args, logger):
    out = open(args.output, 'w')
    chrom_mods = build_transform(args.mod, logger, args.reverse)
    keys = 'chrom chromStart chromEnd name score strand'
    curr_chrom = ""
    with open(args.input, 'r') as input_:
        for line in input_:
            if not line.startswith('#'):
                gene = line.rstrip()
                data = {k:g for k,g in zip(keys.split(), gene.split('\t'))}
                if data.get('chrom') != curr_chrom:
                    curr_chrom = data.get('chrom')
                    positions, deltas = zip(*chrom_mods.get(curr_chrom))
                start_delta = find_delta(positions, 
                                         deltas, 
                                         int(data.get('chromStart')))
                end_delta = find_delta(positions, 
                                       deltas,
                                       int(data.get('chromEnd')))
                data['chromStart'] = int(data.get('chromStart')) + start_delta
                data['chromEnd'] = int(data.get('chromEnd')) + end_delta
                out.write('%s\n' % \
                              '\t'.join([str(data.get(k)) for k in keys.split()]))

def smrna_bed(args, logger):

    out = open(args.output, 'w')
    chrom_mods = build_transform(args.mod, logger, args.reverse)
    keys = 'chrom chromStart chromEnd name score strand'
    curr_chrom = ""
    with open(args.input, 'r') as input_:
        for line in input_:
            gene = line.rstrip()
            data = {k:g for k,g in zip(keys.split(), gene.split('\t'))}
            if data.get('chrom') != curr_chrom:
                curr_chrom = data.get('chrom')
                positions, deltas = zip(*chrom_mods.get(curr_chrom))
            start_delta = find_delta(positions, 
                                     deltas, 
                                     int(data.get('chromStart')))
            end_delta = find_delta(positions, 
                                   deltas,
                                   int(data.get('chromEnd')))
            data['chromStart'] = int(data.get('chromStart')) + start_delta
            data['chromEnd'] = int(data.get('chromEnd')) + end_delta
            out.write('%s\n' % \
                      '\t'.join([str(data.get(k)) for k in keys.split()]))

def bam(args, logger):
    chrom_mods = build_transform(args.mod, logger, args.reverse)
    keys = ''
    curr_chrom = ""

    input_ = AlignmentFile(args.input, 'rb')
    output = AlignmentFile(args.output, 'wb', header=input_.header)

    for line in input_:
        print line.seq, line.cigarstring
