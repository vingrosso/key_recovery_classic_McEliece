from argparse import ArgumentParser
import coloredlogs
from datetime import datetime as dt
from datetime import timedelta as td
import galois
from itertools import combinations
import logging
import msgpack
import numpy as np
from pathlib import Path
from tqdm import tqdm

logger = logging.getLogger(__name__)
coloredlogs.install(logger=logger,
                    level="INFO",
                    datefmt="%H:%M:%S",
                    fmt="[%(asctime)s] %(message)s")

HW = np.array([bin(i).count("1") for i in range(2**13)], dtype="uint32")

Classic_McEliece_polys = {9:  [1, 0, 0, 0, 0, 0, 0, 0, 1, 1],
                          10: [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1],
                          12: [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1],
                          13: [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1]}

def F2m_to_F2(M):
    t, n = np.shape(M)
    return np.reshape(np.swapaxes(M.vector(), 1, 2), (-1, n))


def F2_to_F2m(M, m):
    GF2m = galois.GF(2**m, irreducible_poly=Classic_McEliece_polys[m])
    mt, n = np.shape(M)
    t = mt // m
    return GF2m.Vector(np.swapaxes(np.reshape(M, (t, m, -1)), 2, 1))


def hash_val_to_key(val):
    return bytes([i for i in val])


def keygen(n, t, m):
    GF2m = galois.GF(2**m, irreducible_poly=Classic_McEliece_polys[m])
    while True:
        L = GF2m(np.random.choice(np.arange(2**m), size=n, replace=False))
        while True:
            g = galois.Poly.Random(t, field=GF2m)
            if g.is_irreducible():
                break
        g = g // g.coeffs[0] # monic polynomial
        X = (L[:, np.newaxis]**np.arange(t)).T
        Y = GF2m(np.diag(GF2m.Ones(n) / g(L)))
        H_F2m = X @ Y
        H_F2 = F2m_to_F2(H_F2m)

        Id = np.identity(m * t)
        if np.linalg.matrix_rank(H_F2) == m*t:
            H_pub = H_F2.row_reduce()
            col_Hpub, col_Id = 0, 0
            sub_Id, sub_others = [], []
            while col_Id != m*t:
                if np.array_equal(H_pub[:, col_Hpub], Id[:, col_Id]):
                    sub_Id.append(col_Hpub)
                    col_Id += 1
                else:
                    sub_others.append(col_Hpub)
                col_Hpub += 1
            H_pub[:, np.arange(col_Hpub)] = H_pub[:, sub_Id + sub_others]
            H_F2m[:, np.arange(col_Hpub)] = H_F2m[:, sub_Id + sub_others]
            L[np.arange(col_Hpub)] = L[sub_Id + sub_others]
            assert np.array_equal(H_pub[:, :m*t], galois.GF2.Identity(m*t))
            break
        else:
            logger.info("  Running Keygen again until Hpub = [I|T]")
    return g, L, H_pub, H_F2m


def simulate_HW_leakage(t, m, g, L, accuracy):
    def simulate_perfect_HW_leakage(t, m, g, L):
        GF2m = galois.GF(2**m, irreducible_poly=Classic_McEliece_polys[m])
        leak = np.divide(L[:, np.newaxis]**np.arange(2*t),
                         (g**2)(L)[:, np.newaxis])
        return HW[leak]

    def add_noise(HW_traces, accuracy):
        noise = np.random.choice([-1, 0, 1],
                                 size=np.shape(HW_traces),
                                 p=[(1 - accuracy) / 2,
                                    accuracy,
                                    (1 - accuracy) / 2])
        return np.abs(HW_traces + noise)  # Avoid -1s

    assert (0 < accuracy <= 1)
    perfect_HW_leakage = simulate_perfect_HW_leakage(t, m, g, L)
    if accuracy < 1:
        return add_noise(perfect_HW_leakage, accuracy)
    return perfect_HW_leakage


def load_or_precompute_and_store_hash_table(m, n_powers):
    
    GF2m = galois.GF(2**m, irreducible_poly=Classic_McEliece_polys[m])

    def precompute_hash_table(m):
        hash_table = dict()
        alphas = GF2m(np.arange(1, 2**m, dtype=np.uint16))
        betas = GF2m(np.arange(1, 2**m, dtype=np.uint16))
        alphas_pow_i = alphas[:, np.newaxis]**np.arange(n_powers)
        HW_beta_alphas_pow_i = HW[betas * alphas_pow_i[:, :, np.newaxis]]
        for alpha, HW_seq_plane in tqdm(zip(alphas, HW_beta_alphas_pow_i), total=2**m - 1, desc="Precomputing the hash table"):
            for beta, HW_seq in zip(betas, HW_seq_plane.T):
                key = hash_val_to_key(HW_seq)
                try:
                    multiplicity = int.from_bytes(hash_table[key][-2:], "little")
                    hash_table[key] = bytes(np.array([0, 0, multiplicity + 1], dtype="uint16"))
                except KeyError:
                    hash_table[key] = bytes(np.array([alpha, beta, 1], dtype="uint16"))
        keys_to_remove = [k for k, v in hash_table.items() if int.from_bytes(v[-2:], "little") > 1]
        for key_to_remove in keys_to_remove:
            del hash_table[key_to_remove]
        return hash_table

    hash_table_filename = f"F2^{m}_n_pow_{n_powers}.msgpack"
    if Path(hash_table_filename).is_file():
        logger.info(f"Loading precomputed hash table for F2^{m} and {n_powers=}")
        with open(hash_table_filename, "rb") as hash_table_file:
            h = msgpack.unpackb(hash_table_file.read())
        logger.info(f"Loading done")
    else:
        logger.info(f"Precomputing hash table for F2^{m} and {n_powers=}")
        h = precompute_hash_table(m)
        with open(hash_table_filename, "wb") as hash_table_file:
            hash_table_file.write(msgpack.packb(h))
        logger.info(f"Stored precomputed hash table for F2^{m} and {n_powers=}")
    return h


def attack(H_pub, HW_hash_table, HWs_SCA, t, m, n_powers, delta, tests):

    GF2m = galois.GF(2**m, irreducible_poly=Classic_McEliece_polys[m])

    def match_HW_lists(HW_hash_table, HWs_SCA, t, m, n_powers, delta, tests):
        matches = []
        for i, HWs_seq in enumerate(HWs_SCA):
            key = hash_val_to_key(HWs_seq)
            try:
                raw_bytes = HW_hash_table[key]
                a = int.from_bytes(raw_bytes[0:2], "little")
                gm2a = int.from_bytes(raw_bytes[2:4], "little")
                multiplicity = int.from_bytes(raw_bytes[4:6], "little")
                if multiplicity == 1:
                    matches.append((a, gm2a, i))
                if len(matches) == m * t + delta:
                    logger.info(f"  HW sequences matching completed after {i} iterations, still {n-i} left")
                    return GF2m(matches).T
            except KeyError:
                logger.debug(f"No match for {HWs_seq}")
                pass
        logger.error(f"  Could not match enough HW lists, only {len(matches)}/{m*t+delta}. Possible fixes are:")
        logger.error(f"    - decrease delta (currently {delta}),")
        logger.error(f"    - get a more accurate side-channel distinguisher.")
        exit()


    def interpolate_g(L_a_gm2a_i, m, t, delta, times, tests):
        L_a, L_gm2a, i = L_a_gm2a_i
        L_ga = (GF2m.Ones(t+1) / L_gm2a)**(2**(m-1))
        times.append(dt.now())
        logger.info("  [Alg.1 l.2] Recovered mt+δ pairs (ɑ, g(ɑ))")
        assert(all([ga == tests["g"](a) for (a, ga) in zip(L_a, L_ga)]))
        return galois.lagrange_poly(L_a, L_ga), times


    def recover_permuted_support(H_pub, g, L_a_gm2a_i, t, m, times, tests):
        L_a = L_a_gm2a_i[0, :]
        L_i = np.array([int(i) for i in L_a_gm2a_i[2, :]])
        logger.info(f"    Searching for an invertible submatrix in Hpub")
        for rep, columns_combination in enumerate(combinations(range(m*t+delta), m*t)):
            if np.linalg.matrix_rank(H_pub[:, L_i[list(columns_combination)]]) == m*t:
                break
        logger.info(f"    Invertible submatrix found in Hpub")
        columns_combination = list(columns_combination)
        L_a_sub = L_a[columns_combination]
        assert np.array_equal(L_a_sub, tests['L'][L_i[list(columns_combination)]])
        X = (L_a_sub[:, np.newaxis]**np.arange(t)).T
        Y = GF2m(np.diag(GF2m.Ones(len(L_a_sub)) / g(L_a_sub)))
        V_F2m = X @ Y
        times.append(dt.now())
        logger.info("  [Alg.1 l.4] Vandermonde matrix V constructed with g and the mt+δ (ɑ, g⁻²(ɑ)) pairs")
        assert np.array_equal(V_F2m, tests["H_F2m"][:, L_i[list(columns_combination)]])

        V_F2 = F2m_to_F2(V_F2m)
        assert np.linalg.matrix_rank(V_F2) == m*t
        assert np.linalg.matrix_rank(H_pub[:, L_i[list(columns_combination)]]) == m*t

        S_inv = np.linalg.inv(H_pub[:, L_i[list(columns_combination)]]) @ H_pub
        times.append(dt.now())
        logger.info("  [Alg.1 l.5] Computed the change-of-basis S using V and H_pub")
        H_priv_F2 = V_F2 @ S_inv
        H_priv_F2m = F2_to_F2m(H_priv_F2, m)
        times.append(dt.now())
        logger.info("  [Alg.1 l.6] H_priv = S⁻¹H_pub recovered ")
        return H_priv_F2m[1, :] / H_priv_F2m[0, :], times

    times = [dt.now()]
    L_a_gm2a_i = match_HW_lists(HW_hash_table, HWs_SCA[:, :n_powers], t, m, n_powers, delta, tests)
    g, times = interpolate_g(L_a_gm2a_i[:, :t+1], m, t, delta, times, tests)
    times.append(dt.now())
    logger.info("  [Alg.1 l.3] Polynomial g recovered with interpolation")
    L, times = recover_permuted_support(H_pub, g, L_a_gm2a_i, t, m, times, tests)
    times.append(dt.now())
    for t0, t1 in zip(times[:-1], times[1:]):
        print((t1 - t0) / td(seconds=1), end=", ")
    print()
    logger.info("  [Alg.1 l.7] Full permuted support L recovered ")
    return g, L


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("n",
                        type=int,
                        help="Number of elements in the support L")
    parser.add_argument("t",
                        type=int,
                        help="Error weight")
    parser.add_argument("m",
                        type=int,
                        help="Field degree: GF(2^m)")
    parser.add_argument("--accuracy",
                        type=float,
                        nargs="?",
                        default=1,
                        help="Side-channel distinguisher accuracy")
    args = parser.parse_args()

    n, t, m = args.n, args.t, args.m
    delta = 10
    n_powers = 22

    logger.info("Keygen starts")
    g, L, H_pub, H_F2m = keygen(n, t, m)
    logger.info("Keygen done")
    HWs_SCA = simulate_HW_leakage(t, m, g, L, args.accuracy)
    logger.info(f"Done simulating side-channel trace of shape (n, 2t) = {np.shape(HWs_SCA)}")
    HW_hash_table = load_or_precompute_and_store_hash_table(m, n_powers)
    tests = {"g": g,
             "L": L,
             "H_F2m": H_F2m}
    logger.info("Actual attack starts !")
    t0 = dt.now()
    g_recovered, L_recovered = attack(H_pub, HW_hash_table, HWs_SCA, t, m, n_powers, delta, tests)
    assert g_recovered == g
    assert np.array_equal(L_recovered, L)
    logger.info(f"Attack successful ! Done in {dt.now() - t0}")
