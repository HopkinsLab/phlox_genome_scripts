#!/usr/bin/env python
#
# Script for calculating FIS but really for calculating windowed
# partitioned variances (a,b,c) values based on Weir & Cockerham 1984
# in a pixy-like (https://pixy.readthedocs.io/en/latest/) manner.
# Necessary becasue pixy doesn't give you the raw a, b, and c values.
#
# Note that 'no_snps' here really means number of sites at which there
# is data, which is different in meaning from 'no_snps' in pixy (number of
# variable sites).

import numpy as np
from datatable import dt, f
import allel, argparse
from multiprocessing import Pool
from itertools import combinations

#####################
def contig_table_from_header(vcf):
    # Load contig information by parsing the vcf header
    hdrs = allel.read_vcf_headers(vcf)

    contig_txt = ""
    for line in hdrs.headers:
        if line.startswith('##contig'):
            line = line.replace('##contig=<ID=', '')
            line = line.replace('length=', '')
            line = line.replace('>', '')
            contig_txt += line

    contig_dat = dt.fread(text=contig_txt, sep=',', header=False)
    contig_dat.names = ['CHROM', 'SIZE']
    return contig_dat

def read_popmap(fn):
    # Read species\tpopulations mapping file
    pop_dat = dt.fread(fn, header=False)
    pop_dat.names = ['SP', 'POP']
    return pop_dat

def pop_pair_iter(vcf, populations):
    # Load samples
    samples = np.array(allel.read_vcf_headers(vcf).samples)
    
    # Load population map
    pop_dat = read_popmap(populations)
    uniq_pops = dt.unique(pop_dat['POP']).to_numpy().flatten()
    
    # Get sample indices of each population
    samp_idx = []
    for p in uniq_pops:
        cur_samps = pop_dat[f.POP == p, 'SP'].to_numpy().flatten()
        samp_idx.append(np.where(np.isin(samples, cur_samps))[0].tolist())

    # Return list of indices of population pairs and list of corresponding population names
    return {'idx': list(combinations(samp_idx, 2)), 'popnames': list(combinations(uniq_pops, 2))}

def cut_pos(chunk, max_pos, window_size):
    # Args: VCF chunk (from allel.read_vcf), contig length, and window size
    # Returns dictionary:
    #   start: window start pos (1-based)
    #   end: window end pos (1-based)
    #   idx: list of np arrays specifying indices of chunk's position array that correspond to each window
    lft = list(range(0, max_pos, window_size))
    rgt = lft[1:] + [max_pos]
    out_dict = {'start': np.array(lft)+1, 'end': np.array(rgt)}
    if chunk is None:
        out_dict['idx'] = [ np.array([], dtype='int64') for i in range(len(lft))]
    else:
        pos_arr = chunk['variants/POS']
        window_idx = [np.where(np.logical_and(pos_arr > x, pos_arr <= y))[0] for x, y in zip(lft, rgt)]
        out_dict['idx'] = window_idx
    return out_dict

def compute_abc_on_window(g, subpops):
    mask = np.logical_and(
        np.any(np.any(g[:, subpops[0], :] >= 0, axis=1), axis=1),
        np.any(np.any(g[:, subpops[1], :] >= 0, axis=1), axis=1)
    )
    no_snps = mask.sum()
    if no_snps > 0:
        a, b, c = allel.weir_cockerham_fst(g[mask, :, :], subpops)
        out_arr = np.array([np.nansum(a), np.nansum(b), np.nansum(c), no_snps])
    else:
        out_arr = np.array([np.nan, np.nan, np.nan, 0])
    return out_arr

def write_header(filename, hdr_cols):
    with open(filename, 'wt') as fh:
        fh.write('\t'.join(hdr_cols) + '\n')
        
def write_output(filename, out_dat):
    if out_dat.shape[0] > 0:
        out_dat.to_pandas().to_csv(filename, sep="\t", mode='a', header=False, index=False, na_rep='NA')

#####################
def main():
    ### Program args
    prog_description = ('Calculates FIS using formulae from Weir and Cockerham 1984.'
                        '')

    parser = argparse.ArgumentParser(
                        prog='FIS.py',
                        description=prog_description,
                        epilog='')

    parser.add_argument('--vcf'        , type=str, required=True)
    parser.add_argument('--populations', type=str, required=True)
    parser.add_argument('--outfile',     type=str, required=True)
    parser.add_argument('--chromosomes', type=str, required=False)
    parser.add_argument('--window_size', type=int, default=10000)
    parser.add_argument('--n_cores',     type=int, default=1)

    args = parser.parse_args()

    
    vcf = args.vcf
    populations = args.populations

    ### Load contig information
    contig_dat = contig_table_from_header(vcf)
    if args.chromosomes is not None:
        # Filter chromosomes if need be.
        keep_chrom = [ x for x in args.chromosomes.split(',') if x != '']
        mask = dt.rowany(contig_dat[:, [f.CHROM == x for x in keep_chrom]])
        contig_dat = contig_dat[mask,:]
        if contig_dat.shape[0] == 0:
            raise ValueError("None of chromosomes listed are in the VCF header")
    print(f'Analyzing {contig_dat.shape[0]} contigs:')

    ### Get population pairs and iterator
    pop_pairs = pop_pair_iter(vcf, populations)

    ### Write header of output
    out_hdr_cols = ['pop1', 'pop2', 'chromosome', 'window_pos_1', 'window_pos_2', 'a', 'b', 'c', 'avg_fis', 'avg_fst', 'avg_fit', 'no_snps']
    write_header(args.outfile, out_hdr_cols)

    ### To avoid reading everything into memory at once, use tabix to grab each chromosome
    with Pool(args.n_cores) as pool_obj:
        for chrom, csize in contig_dat.to_tuples():
            print(f'  Processing contig {chrom} [{csize}] ...')
            chunk = allel.read_vcf(
                input=vcf, 
                fields=['variants/CHROM', 'variants/POS', 'calldata/GT'],
                region=chrom,
                tabix='tabix')
            npos = chunk['variants/POS'].shape[0]
            print(f'    npos = {npos}')

            window_info = cut_pos(chunk, csize, args.window_size)

            npairs = len(pop_pairs['idx'])
            nwindows = len(window_info['idx'])

            # Step through every 50 pairs of populations
            stepsize = min(npairs, 50)
            step = 0
            while step < npairs:
                pair_range = range(step, min(npairs, step + stepsize))
                glist = [chunk['calldata/GT'][x,:,:] for pair_i in pair_range for x in window_info['idx']]
                plist = [pop_pairs['idx'][pair_i]    for pair_i in pair_range for i in range(nwindows)]
                
                out_dat = dt.Frame(np.array(pool_obj.starmap(compute_abc_on_window, zip(glist, plist))))
                
                out_dat.names = ['a', 'b', 'c', 'no_snps']
                out_dat['no_snps'] = dt.int32
                # out_dat['pop1'] = pop_pairs['popnames'][pair_i][0]
                # out_dat['pop2'] = pop_pairs['popnames'][pair_i][1]
                out_dat['pop1'] = np.array([pop_pairs['popnames'][pair_i][0] for pair_i in pair_range for i in range(nwindows)])
                out_dat['pop2'] = np.array([pop_pairs['popnames'][pair_i][1] for pair_i in pair_range for i in range(nwindows)])
                out_dat['chromosome'] = chrom
                # out_dat['window_pos_1'] = window_info['start']
                # out_dat['window_pos_2'] = window_info['end']
                out_dat['window_pos_1'] = np.concatenate([window_info['start'] for pair_i in pair_range])
                out_dat['window_pos_2'] = np.concatenate([window_info['end']   for pair_i in pair_range])
                out_dat[(f.b + f.c) != 0,       dt.update(avg_fis = 1 - (f.c / (f.b + f.c)))]
                out_dat[(f.a + f.b + f.c) != 0, dt.update(avg_fst = f.a / (f.a + f.b + f.c))]
                out_dat[(f.a + f.b + f.c) != 0, dt.update(avg_fit = 1 - f.c / (f.a + f.b + f.c))]
                
                write_output(args.outfile, out_dat[f.no_snps > 0, out_hdr_cols])
                step += stepsize

#####################


if __name__ == '__main__':
    main()

    

