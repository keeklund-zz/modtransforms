def gtf(args):
    out = open(args.output, 'w')

    chrom_mods = build_transform(args.mod)
    keys = 'contig feature source start end score strand frame attributes'
    curr_chrom = ""
    with open(args.input, 'r') as input_:
        for line in input_:
            gene = line.rstrip()
            data = {k:g for k,g in zip(keys.split(), gene.split('\t'))}
            if data.get('contig') != curr_chrom:
                curr_chrom = data.get('contig')
                positions, deltas = zip(*chrom_mods.get(curr_chrom))
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
            except:
                pass

def smrna_gff3(args):
    out = open(args.output, 'w')
    chrom_mods = build_transform(args.mod)
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
#                try:
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
#                except:
#                    pass

def smrna_bed(args):

    out = open(args.output, 'w')
    chrom_mods = build_transform(args.mod)
    keys = 'chrom chromStart chromEnd name score strand'
    curr_chrom = ""
    with open(args.input, 'r') as input_:
        for line in input_:
            gene = line.rstrip()
            data = {k:g for k,g in zip(keys.split(), gene.split('\t'))}
            if data.get('chrom') != curr_chrom:
                curr_chrom = data.get('chrom')
                positions, deltas = zip(*chrom_mods.get(curr_chrom))
            try:
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
            except:
                pass

def file_type_error(args):
    exit("ERROR: '%s' is not a valid input type" % args.type)

