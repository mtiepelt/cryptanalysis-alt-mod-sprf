"""
    Figure 2: distribution of Hamming weights.
"""
import argparse
import sys

from paths import SRC_DIR
sys.path.append(str(SRC_DIR))

import parameters
import distributionHW as dist

# Values used in the paper (--reproduce-paper)
PAPER_NUM_EXPS = 100
PAPER_AN = 4.4
PAPER_Q = 5
PAPER_T = 2**19
PAPER_CONST = 0.9


def run(secpar, num_exps, an, q, T, const, num_cpus=None):
    """
        One run for the given secpar, with the largest al that satisfies the
        success condition.
    """
    # Be aware: SageMath 'right_kernel()' is memory hungry!
    al = parameters.bin_search_largest_l(secpar, an, q, T, T, const,
                                         prec=1000)
    dist.stats(secpar, an, al, q, num_exps, T=T, num_cpus=num_cpus)


def reproduce_paper(num_cpus=None):
    """
        Computations for Figure 2 in the paper, with the paper's values.
    """
    num_exps = PAPER_NUM_EXPS
    an = PAPER_AN
    q = PAPER_Q
    T = PAPER_T
    const = PAPER_CONST

    # Be aware: SageMath 'right_kernel()' is memory hungry!
    dist.stats(20, an,
               parameters.bin_search_largest_l(20, an, q, T, T, const,
                                               prec=1000),
               q, num_exps, T=T, num_cpus=num_cpus)
    dist.stats(64, an,
               parameters.bin_search_largest_l(64, an, q, T, T, const,
                                               prec=1000),
               q, num_exps, T=T, num_cpus=num_cpus)
    dist.stats(128, an,
               parameters.bin_search_largest_l(128, an, q, T, T, const,
                                               prec=1000),
               q, num_exps, T=T, num_cpus=num_cpus)


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n', '--secpar', type=int,
                        help='security parameter '
                             '(required unless --reproduce-paper)')
    parser.add_argument('--reproduce-paper', action='store_true',
                        help='run the paper\'s computations with the paper\'s '
                             'values (secpar = 20, 64, 128, num_exps = 100, '
                             'an = 4.4, q = 5, T = 2**19, const = 0.9); '
                             'all other options except --num-cpus are ignored')
    parser.add_argument('--num_exps', type=int, default=100,
                        help='number of experiments')
    parser.add_argument('--an', type=float, default=4.4,
                        help='n = ceil(an * secpar)')
    parser.add_argument('--q', type=int, default=5,
                        help='d = n - l + q')
    parser.add_argument('--log2T', type=int, default=None,
                        help='T = 2**log2T; None means secpar - 1')
    parser.add_argument('--num-cpus', type=int, default=None,
                        help='number of CPUs; None uses all but two')
    parser.add_argument('--const', type=float, default=0.9,
                        help='target success probability in the success '
                             'condition (Eq. 3)')
    args = parser.parse_args()

    if args.reproduce_paper and args.secpar is not None:
        parser.error('--secpar cannot be combined with --reproduce-paper')
    if not args.reproduce_paper and args.secpar is None:
        parser.error('--secpar is required unless --reproduce-paper is given')
    return args


if __name__ == '__main__':
    args = parse_args()

    if args.reproduce_paper:
        reproduce_paper(args.num_cpus)
    else:
        log2T = args.secpar - 1 if args.log2T is None else args.log2T
        run(args.secpar, args.num_exps, args.an, args.q, 2**log2T, args.const,
            args.num_cpus)
