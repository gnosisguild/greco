import os
import json
from bfv.crt import CRTModuli
from bfv.bfv import BFVCrt
from bfv.discrete_gauss import DiscreteGaussian
from bfv.polynomial import Polynomial, poly_div


def main(args):
    n = args.n
    qis = json.loads(args.qis)
    t = args.t

    crt_moduli = CRTModuli(qis)
    sigma = 3.2
    discrete_gaussian = DiscreteGaussian(sigma)
    bfv_crt = BFVCrt(crt_moduli, n, t, discrete_gaussian)

    s = bfv_crt.SecretKeyGen()
    e = bfv_crt.bfv_q.rlwe.SampleFromErrorDistribution()
    ais = [bfv_crt.bfv_qis[i].rlwe.Rq.sample_polynomial() for i in range(len(crt_moduli.qis))]

    pub_key = bfv_crt.PublicKeyGen(s, e, ais)

    m = bfv_crt.bfv_q.rlwe.Rt.sample_polynomial()
    e0 = bfv_crt.bfv_q.rlwe.SampleFromErrorDistribution()
    e1 = bfv_crt.bfv_q.rlwe.SampleFromErrorDistribution()
    u = bfv_crt.bfv_q.rlwe.SampleFromTernaryDistribution()

    ciphertext = bfv_crt.PubKeyEncrypt(pub_key, m, e0, e1, u)

    k1 = m.scalar_mul(crt_moduli.q)
    k1.reduce_coefficients_by_modulus(t)

    p = 21888242871839275222246405745257275088548364400416034343698204186575808495617

    r2is = []
    r1is = []
    ct0is = []
    ct1is = []
    pk0is = []
    pk1is = []
    p1is = []
    p2is = []

    pk_array = [pk[0].coefficients for pk in pub_key]
    pk1_array = [pk[1].coefficients for pk in pub_key]

    cyclo = [1] + [0] * (n - 1) + [1]
    cyclo = Polynomial(cyclo)

    for i, cti in enumerate(ciphertext):
        ct0i = cti[0]
        ct1i = cti[1]

        k0i = pow(-t, -1, qis[i])

        pk0 = Polynomial(pk_array[i])
        pk1 = Polynomial(pk1_array[i])

        ct0i_hat = pk0 * u + e0 + k1.scalar_mul(k0i)
        num = ct0i + ct0i_hat.scalar_mul(-1)
        num.reduce_coefficients_by_modulus(qis[i])
        quotient, rem = poly_div(num.coefficients, cyclo.coefficients)
        r2i = Polynomial(quotient)

        num = ct0i + ct0i_hat.scalar_mul(-1) + (r2i * cyclo).scalar_mul(-1)
        quotient, rem = poly_div(num.coefficients, [qis[i]])
        r1i = Polynomial(quotient)

        ct1i_hat = pk1 * u + e1

        num = ct1i + ct1i_hat.scalar_mul(-1)
        num.reduce_coefficients_by_modulus(qis[i])
        quotient, rem = poly_div(num.coefficients, cyclo.coefficients)
        p2i = Polynomial(quotient)

        num = ct1i + ct1i_hat.scalar_mul(-1) + (p2i * cyclo).scalar_mul(-1)
        quotient, rem = poly_div(num.coefficients, [qis[i]])
        p1i = Polynomial(quotient)

        pk1is.append(pk1)
        p2is.append(p2i)
        p1is.append(p1i)
        ct1is.append(ct1i)
        r2is.append(r2i)
        r1is.append(r1i)
        ct0is.append(ct0i)
        pk0is.append(pk0)


    pk0i_assigned = [Polynomial([(coef + p) % p for coef in pk0i.coefficients]) for pk0i in pk0is]
    pk1i_assigned = [Polynomial([(coef + p) % p for coef in pk1i.coefficients]) for pk1i in pk1is]

    e0_assigned = Polynomial([(coef + p) % p for coef in e0.coefficients])
    e1_assigned = Polynomial([(coef + p) % p for coef in e1.coefficients])
    k1_assigned = Polynomial([(coef + p) % p for coef in k1.coefficients])
    u_assigned = Polynomial([(coef + p) % p for coef in u.coefficients])

    r1is_assigned = [Polynomial([(coef + p) % p for coef in r1i.coefficients]) for r1i in r1is]
    r2is_assigned = [Polynomial([(coef + p) % p for coef in r2i.coefficients]) for r2i in r2is]
    p1is_assigned = [Polynomial([(coef + p) % p for coef in p1i.coefficients]) for p1i in p1is]
    p2is_assigned = [Polynomial([(coef + p) % p for coef in p2i.coefficients]) for p2i in p2is]
    ct0is_in_p = [Polynomial([(coef + p) % p for coef in ct0i.coefficients]) for ct0i in ct0is]
    ct1is_in_p = [Polynomial([(coef + p) % p for coef in ct1i.coefficients]) for ct1i in ct1is]

    json_input = {
        "pk0_qi": [[str(coef) for coef in pk0i.coefficients] for pk0i in pk0i_assigned],
        "pk1_qi": [[str(coef) for coef in pk1i.coefficients] for pk1i in pk1i_assigned],
        "u": [str(coef) for coef in u_assigned.coefficients],
        "e0": [str(coef) for coef in e0_assigned.coefficients],
        "e1": [str(coef) for coef in e1_assigned.coefficients],
        "k1": [str(coef) for coef in k1_assigned.coefficients],
        "r2is": [[str(coef) for coef in r2i.coefficients] for r2i in r2is_assigned],
        "r1is": [[str(coef) for coef in r1i.coefficients] for r1i in r1is_assigned],
        "p2is": [[str(coef) for coef in p2i.coefficients] for p2i in p2is_assigned],
        "p1is": [[str(coef) for coef in p1i.coefficients] for p1i in p1is_assigned],
        "ct0is": [[str(coef) for coef in ct0i_in_p.coefficients] for ct0i_in_p in ct0is_in_p],
        "ct1is": [[str(coef) for coef in ct0i_in_p.coefficients] for ct0i_in_p in ct1is_in_p],
    }

    qis_bitsize = max(qis).bit_length()
    qis_len = len(qis)
    filename = f"pk_enc_{args.n}_{qis_len}x{qis_bitsize}_{args.t}.json"
    output_path = os.path.join("src", "data", "pk_enc_data", filename)

    with open(output_path, 'w') as f:
        json.dump(json_input, f)

    json_input_zeroes = {
        "pk0_qi": [["0" for _ in pk0i.coefficients] for pk0i in pk0i_assigned],
        "pk1_qi": [["0" for _ in pk1i.coefficients] for pk1i in pk1i_assigned],
        "u": ["0" for _ in u_assigned.coefficients],
        "e0": ["0" for _ in e0_assigned.coefficients],
        "e1": ["0" for _ in e1_assigned.coefficients],
        "k1": ["0" for _ in k1_assigned.coefficients],
        "r2is": [["0" for _ in r2i.coefficients] for r2i in r2is_assigned],
        "r1is": [["0" for _ in r1i.coefficients] for r1i in r1is_assigned],
        "p2is": [["0" for _ in p2i.coefficients] for p2i in p2is_assigned],
        "p1is": [["0" for _ in p1i.coefficients] for p1i in p1is_assigned],
        "ct0is": [["0" for _ in ct0i_in_p.coefficients] for ct0i_in_p in ct0is_in_p],
        "ct1is": [["0" for _ in ct1i_in_p.coefficients] for ct1i_in_p in ct1is_in_p],
    }

    output_path = os.path.join("src", "data", "pk_enc_data", f"pk_enc_{args.n}_{qis_len}x{qis_bitsize}_{args.t}_zeroes.json")

    with open(output_path, 'w') as f:
        json.dump(json_input_zeroes, f)



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate rust constants and json inputs for BFV zk proof of secret key encryption circuit"
    )
    parser.add_argument(
        "-n", type=int, required=True, help="Degree of f(x), must be a power of 2."
    )
    parser.add_argument(
        "-qis", type=str, required=True, help="List of qis such that qis[i] is the modulus of the i-th CRT basis of the modulus q of the ciphertext space."
    )
    parser.add_argument(
        "-t", type=int, required=True, help="Modulus t of the plaintext space."
    )

    args = parser.parse_args()
    main(args)
