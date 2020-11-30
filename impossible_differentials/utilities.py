def bits(n, size):
    output = [0] * size
    for i in range(size):
        output[i] = (n >> i) & 1
    return output


def log2(n):
    return n.bit_length() - 1


def ddt_rows(ddt, in_size, out_size):
    out = {}
    for a in range(1 << in_size):
        out[a] = set()
        for b in range(1 << out_size):
            if ddt[a, b] != 0:
                out[a].add(b)
    return out


def ddt_cols(ddt, in_size, out_size):
    out = {}
    for b in range(1 << out_size):
        out[b] = set()
        for a in range(1 << in_size):
            if ddt[a, b] != 0:
                out[b].add(a)
    return out


# def product(my_list, repeat=1):
#     # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
#     # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
#     pools = [tuple(pool) for pool in my_list] * repeat
#     result = [[]]
#     for pool in pools:
#         result = [x + [y] for x in result for y in pool]
#     for prod in result:
#         yield tuple(prod)


def hwt(x):
    return bin(x).count("1")
