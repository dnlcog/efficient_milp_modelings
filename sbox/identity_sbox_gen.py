import pickle
import sys


def generate_identity_sbox(size):
    ineg_set = set()
    for i in range(size):
        ineg = [0] * (2 * size + 1)
        ineg[i] = 1
        ineg[i + size] = -1
        ineg_set.add(tuple(ineg))
        ineg = [0] * (2 * size + 1)
        ineg[i] = -1
        ineg[i + size] = 1
        ineg_set.add(tuple(ineg))

    ddt = {}
    ddt[0, 0] = 1
    for a in range(1, 1 << size):
        for b in range(1 << size):
            ddt[a, b] = 0 if a != b else 1
        ddt[0, a] = 0

    with open("identity_sbox_{}.pkl".format(size), "wb") as f:
        pickle.dump((size, size, ddt, ineg_set), f, 3)


if __name__ == "__main__":
    try:
        size = int(sys.argv[1])
    except IndexError:
        raise SystemExit("Usage: {} ".format(sys.argv[0]) + "<size>")

    generate_identity_sbox(size)
