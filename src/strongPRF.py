from sage.all import GF, matrix, vector

import pcommon


def vec_F3_to_F2(vec, sparse=False):
    # Modulus transformation
    if sparse:
        # Sparse: Maps {0,2} --> {0}, {1} --> {1}
        conv_map = {0: 0, 1: 1, 2: 0}
    else:
        # Non-Sparse: Maps {0} --> {0}, {1,-1} --> {1}
        conv_map = {0: 0, 1: 1, 2: 1}

    return [conv_map[coef] for coef in vec]


def G_x_F3(G, x):
    """
        G * (x||1) over F3.
    """
    return G * vector(GF(3), list(x) + [1])


def G_x_F2(G, x, sparse):
    """
        G * (x||1) over F3, mapped to F2.
    """
    return vec_F3_to_F2(G_x_F3(G, x), sparse)


class simple_sPRF:
    """
        (2,3) Strong PRF following
        [Def. 9, https://eprint.iacr.org/2024/582.pdf]

        Note: here q is the 2nd modulus, not a parameter of the attack.
    """

    def __init__(self, secpar, d, n, m, t, p=2, q=3, sparse=False, seed=None):
        """
            seed: optional, makes G, A, B and the key reproducible.
        """
        rng = pcommon.new_rng(seed)

        self.d = d
        self.n = n
        self.m = m
        self.t = t

        self.secpar = secpar
        self.sparse = sparse

        F3 = GF(q)
        F2 = GF(p)

        self.G = matrix(F3, rng.integers(0, q, size=(n, d + 1)))
        self.A = matrix(F2, rng.integers(0, p, size=(m, n)))
        self.B = matrix(F3, rng.integers(0, q, size=(t, m)))
        self.key = vector(F2, rng.integers(0, p, size=(n)))

    def eval_sPRF(self, x):
        z = G_x_F2(self.G, x, self.sparse)

        # wPRF(z)
        x_k = [ki * zi for ki, zi in zip(self.key, z)]
        A_x_k = self.A * vector(GF(2), x_k)
        B_A_x_k = self.B * vector(GF(3), [i for i in A_x_k])
        return B_A_x_k
