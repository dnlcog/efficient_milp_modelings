def bits(n, size):
    output = [0] * size
    for i in range(size):
        output[i] = (n >> i) & 1
    return output


def log2(n):
    return n.bit_length() - 1


def scalar_prod(a, b):
    res = 0
    for x, y in zip(a, b):
        res += x * y
    return res


def point_kept(point, ineg):
    alpha = scalar_prod(bits(point, len(ineg) - 1), ineg[:-1]) + ineg[-1]
    return alpha >= 0


def vertex_is_on_face(point, ineg):
    alpha = scalar_prod(bits(point, len(ineg) - 1), ineg[:-1]) + ineg[-1]
    return alpha == 0


def xor_cons(model, input_list, output):
    # output represents a bit and input_list
    # several bits
    n = len(input_list)

    for i in range(1 << n):
        bit_list = bits(i, n)
        constant = -1 + sum(bit_list)
        constraint = qs(
            input_list[j] if bit_list[j] == 0 else -input_list[j] for j in range(n)
        )
        xor = 0
        for j in range(n):
            xor ^= bit_list[j]

        if xor == 0:
            constraint -= output
            constant += 1
        else:
            constraint += output


def nor_cons(model, input_list, output):
    n = len(input_list)

    for i in range(1 << n):
        bit_list = bits(i, n)
        constant = -1 + sum(bit_list)
        constraint = qs(
            input_list[j] if bit_list[j] == 0 else -input_list[j] for j in range(n)
        )
        res = 1
        for j in range(n):
            res &= 1 - bit_list[j]

        if res == 0:
            constraint -= output
            constant += 1
        else:
            constraint += output
