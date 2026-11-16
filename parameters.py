import math

from paths import PARA_DIR


def get_l(secpar, al):
    """
        l = ceil(al * secpar).
        Round first, since floating-point products can be slightly off
        (e.g. 0.07 * 100 = 7.000000000000001).
    """
    return math.ceil(round(al * secpar, 9))


def sPRF_parameters(secpar, q, an, al):
    # Round first, since floating-point products can be slightly off
    # (e.g. 0.07 * 100 = 7.000000000000001)
    n = math.ceil(round(an * secpar, 9))
    l = get_l(secpar, al)
    d = n - l + q
    m = math.ceil(7.06 * secpar)
    t = math.ceil(2 / math.log(3, 2) * secpar)
    return n, l, d, m, t


def log3_approxN_closed(N, R, l):
    """
        Expected Hamming weight X_R of the R-th smallest of N Hamming weights
        of uniform vectors in F_3^l (approximation from normal distribution).
    """
    inverfc_bound = - math.sqrt(2 * math.log(N / R)
                                - math.log(4 * math.pi)
                                - math.log(math.log(N / (2 * R))))
    X_R = (2 * l / 3) + math.sqrt(2 * l) / 3 * inverfc_bound - 1 / 2
    return X_R


def get_expected_pairs(d):
    """
        E[2^hw_0(y)]
    """
    return (4 / 3)**d


def get_N(T1):
    return int(T1)


def get_R(T2, d):
    return T2 / get_expected_pairs(d)


def attack_budget(T1, T2, d, l):
    """
        Number of kernel elements N, number of kept kernel elements R
        (rounded up), and X_R: the expected Hamming weight of the R-th
        smallest element (see log3_approxN_closed).
    """
    N = get_N(T1)
    R_int = math.ceil(get_R(T2, d))
    log3_ZR = log3_approxN_closed(N, R_int, l)
    return N, R_int, log3_ZR


def satisfies_success_condition(T1, T2, d, l, const):
    """
        Success condition (Eq. 3)
    """
    N, R_int, log3_ZR = attack_budget(T1, T2, d, l)
    const_factor = -math.log(1 - const)
    return math.log(const_factor, 3) + log3_ZR <= math.log(T2, 3)


def bin_search_largest_l(secpar, an, q, T1, T2, const, prec):
    """
        Binary search for the largest al (in steps of 1/prec) that satisfies
        the success condition (Eq. 3)
    """
    best_al = 0

    # Search over al = v / prec for integers lo <= v <= hi
    lo = 0
    hi = round(an * prec)

    while lo <= hi:
        mid = (lo + hi) // 2
        al = mid / prec

        n, l, d, m, t = sPRF_parameters(secpar, q, an=an, al=al)

        if satisfies_success_condition(T1, T2, d, l, const):
            best_al = al
            lo = mid + 1
        else:
            hi = mid - 1
    return best_al


def bin_search_smallest_logT(secpar, an, q, al, const, prec=10):
    """
        Binary search for the smallest log_2 T (in steps of 1/prec) that
        satisfies the success condition (Eq. 3)
    """
    best_logT = 0

    # Search over log_2 T = v / prec for integers lo <= v <= hi
    lo = 0
    hi = secpar * prec

    n, l, d, m, t = sPRF_parameters(secpar, q, an=an, al=al)

    while lo <= hi:
        mid = (lo + hi) // 2
        T = 2**(mid / prec)

        if satisfies_success_condition(T, T, d, l, const):
            best_logT = mid / prec
            hi = mid - 1
        else:
            lo = mid + 1

    return best_logT


def bounds_l(an, q, T, range_secpar, const):
    """
        Compute beta as a function of T.
        T = None: T1 = T2 = 2**(secpar - 1).
    """
    secpar_to_max_al = {}
    for secpar in range_secpar:
        T_secpar = 2**(secpar - 1) if T is None else T

        max_al = bin_search_largest_l(secpar, an, q, T_secpar, T_secpar,
                                      const, prec=1000)
        secpar_to_max_al[secpar] = max_al

    return secpar_to_max_al


def bounds_T_for_l(secpar, an, q, min_al, const, prec=10, prec_al=1000,
                   step=None):
    """
        Compute T as a function of beta
    """
    al_to_minlogT = {}
    max_al = bin_search_largest_l(secpar, an, q, 2**(secpar - 1),
                                  2**(secpar - 1), const, prec=prec_al)

    if step is None:
        step = int((max_al - min_al) * prec_al // 8)

    for v in range(round(max_al * prec_al), round(min_al * prec_al), -step):
        al = v / prec_al
        log_T2 = bin_search_smallest_logT(secpar, an, q, al, const, prec=prec)
        al_to_minlogT[al] = log_T2  # log_2 T_2, with T_1 = T_2

    return al_to_minlogT


def l_values(secpar_to_al):
    """
        l = ceil(al * secpar) for every (secpar, al), as used in the attack
    """
    return [get_l(secpar, al) for secpar, al in secpar_to_al.items()]


def write_bound(secpar_to_al, fn, l_row=True):
    """
        Export CSV to read into TeX/PGF
    """
    with open(fn, 'w') as f:
        f.write(','.join([str(x) for x in secpar_to_al.keys()]) + '\n')
        f.write(','.join([str(x) for x in secpar_to_al.values()]) + '\n')
        if l_row:
            f.write(','.join([str(x) for x in l_values(secpar_to_al)]))


def print_table(header, rows, width=5):
    """
        Print a table: rows is a list of (label, values).
        A separator line follows the first row.
    """
    label_width = max(len(label) for label, _ in rows)

    print("-" * 50)
    print(f" {header}")
    print("-" * 50)
    for i, (label, values) in enumerate(rows):
        print(label.ljust(label_width) + " : "
              + " ".join(str(x).rjust(width) for x in values))
        if i == 0:
            print("-" * 5)
    print("-" * 50)


if __name__ == '__main__':
    an = 4.4
    const = 0.9

    # Table 4
    # T as function of beta
    secpar = 128
    q = 0
    al_to_minlogT = bounds_T_for_l(secpar, an, q, min_al=0.1, const=const,
                                   prec=10)
    write_bound(al_to_minlogT,
                f'{PARA_DIR}/al-T_bounds_s-{secpar}_an-{an}_q-{q}.tab',
                l_row=False)
    print_table(f"Table 1: secpar = {secpar}",
                [("beta", al_to_minlogT.keys()),
                 ("log_2 T_2", al_to_minlogT.values())])

    # Table 5
    # beta as function of T
    q = 0
    secpar_to_max_al = bounds_l(an, q, None, list(range(32, 256 + 1, 32)),
                                const)
    print()
    write_bound(secpar_to_max_al,
                f'{PARA_DIR}/T-secpar_bounds_an-{an}_q-{q}.tab')
    print_table(f"Table 5: T1 + T2 = T <= 2**secpar",
                [("secpar", secpar_to_max_al.keys()),
                 ("al", secpar_to_max_al.values()),
                 ("l", l_values(secpar_to_max_al))])

    # Table 6
    q = 5
    log2T = 19
    secpar_to_max_al_fixed = bounds_l(an, q, 2**log2T,
                                      [20] + list(range(32, 256 + 1, 32)),
                                      const)
    write_bound(secpar_to_max_al_fixed,
                f'{PARA_DIR}/T-{log2T}_bounds_an-{an}_q-{q}.tab')
    print_table(f"Table 6: T1 = T2 = 2**{log2T}, q={q}",
                [("secpar", secpar_to_max_al_fixed.keys()),
                 ("al", secpar_to_max_al_fixed.values()),
                 ("l", l_values(secpar_to_max_al_fixed))])
