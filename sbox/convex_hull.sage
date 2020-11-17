from utilities import *
from sage.crypto.sboxes import SBox
from sage.crypto.sboxes import sboxes
import itertools
import pickle
import time
import argparse


def sbox(sbox_name):
    """
    Returns the SBox object from its name.
    """
    if sbox_name == "Keccak":
        return SBox(
            [
                0,
                9,
                18,
                11,
                5,
                12,
                22,
                15,
                10,
                3,
                24,
                1,
                13,
                4,
                30,
                7,
                20,
                21,
                6,
                23,
                17,
                16,
                2,
                19,
                26,
                27,
                8,
                25,
                29,
                28,
                14,
                31,
            ]
        )
    elif sbox_name == "Pyj3":
        return SBox([1, 3, 6, 5, 2, 4, 7, 0])
    elif sbox_name == "Pyj4":
        return SBox([2, 0xD, 3, 9, 7, 0xB, 0xA, 6, 0xE, 0, 0xF, 4, 8, 5, 1, 0xC])
    else:
        return sboxes[sbox_name]


def sets_from_ddt(ddt):
    """
    Given a DDT, returns the sets of possible and impossible points.
    """
    pos_trans = set([])
    imp_trans = set([])
    n, m = ddt.dimensions()
    for in_diff in range(n):
        for out_diff in range(m):
            point = int(in_diff + (out_diff << log2(n)))
            if ddt[in_diff, out_diff] == 0:
                imp_trans.add(point)
            else:
                pos_trans.add(point)
    return (pos_trans, imp_trans)


def inequalities(sbox, nb_faces=2):
    """
    Given an SBox object, computes big set of inequalities
    with faces from the convex hull of possible points and
    the additions of at most nb_faces of them.
    """
    ddt = sbox.difference_distribution_table()
    in_size = sbox.input_size()
    out_size = sbox.output_size()
    n = in_size + out_size

    pos_trans, imp_trans = sets_from_ddt(ddt)
    point_set = imp_trans.copy()

    print("Building polyhedron")
    pos_polyhedron = Polyhedron(vertices=[bits(point, n) for point in pos_trans])

    # We collect faces of the convex hull as pure python objects.
    ineqs = [[int(x) for x in ineg] for ineg in pos_polyhedron.inequalities()]

    # Our personal way of writing an inequality puts the constant
    # at the end of the tuple.
    ineqs = [tuple(ineg[1:] + [ineg[0]]) for ineg in ineqs]

    eqs = pos_polyhedron.equations()
    assert len(eqs) == 0

    # We keep the polyhedron faces.
    ineg_set = set([x for x in ineqs])

    print("Computing discarded points for faces.")
    ineg_to_points = {}
    for ineg in ineg_set:
        ineg_to_points[ineg] = set([])
        for point in point_set:
            if not (point_kept(point, ineg)):
                ineg_to_points[ineg].add(point)

    print("Currently {} inequalities. Removing inclusions...".format(len(ineg_set)))
    for ineg1 in ineg_set:
        if ineg1 in ineg_to_points:
            for ineg2 in ineg_set:
                if (ineg1 != ineg2) and (ineg2 in ineg_to_points):
                    if ineg_to_points[ineg2] <= ineg_to_points[ineg1]:
                        del ineg_to_points[ineg2]

    ineg_set = set(ineg_to_points.keys())

    ineqs = ineg_set.copy()

    print("Now {} inequalities.".format(len(ineg_set)))

    if nb_faces >= 2:

        # We add new equations obtained by summing up to nb_faces faces.
        count = 0
        for center_point in pos_trans:
            count += 1
            print("center point: {}".format(count))

            faces = [face for face in ineqs if vertex_is_on_face(center_point, face)]
            print("  # of faces: {}".format(len(faces)))

            local_inegs = set([x for x in faces])
            add_counter = 0
            total_counter = 0

            for face_indices in itertools.product(range(len(faces)), repeat=nb_faces):
                total_counter += 1
                new_ineq = [int(0) for i in range(n + 1)]
                for face_index in face_indices:
                    face = faces[face_index]
                    for i in range(n + 1):
                        new_ineq[i] += face[i]
                new_ineg = tuple(new_ineq)

                if new_ineg not in ineg_set:

                    new_elim_points = set([])
                    for point in point_set:
                        if not (point_kept(point, new_ineg)):
                            new_elim_points.add(point)

                    to_add = True
                    for ineg in local_inegs:
                        if new_elim_points <= ineg_to_points[ineg]:
                            to_add = False
                            break

                    if to_add:
                        add_counter += 1
                        local_inegs.add(new_ineg)
                        ineg_set.add(new_ineg)
                        ineg_to_points[new_ineg] = new_elim_points

            print("    added {} inequalities.".format(add_counter))

        print("Currently {} inequalities. Removing inclusions...".format(len(ineg_set)))
        for ineg1 in ineg_set:
            if ineg1 in ineg_to_points:
                for ineg2 in ineg_set:
                    if (ineg1 != ineg2) and (ineg2 in ineg_to_points):
                        if ineg_to_points[ineg2] <= ineg_to_points[ineg1]:
                            del ineg_to_points[ineg2]

        ineg_set = set(ineg_to_points.keys())

    print("Building point_to_inegs then checking.")
    point_to_inegs = {}
    for point in point_set:
        point_to_inegs[point] = set([])
        for ineg in ineg_set:
            if point in ineg_to_points[ineg]:
                point_to_inegs[point].add(ineg)
        assert len(point_to_inegs[point]) >= 1

    for point in pos_trans:
        for ineg in ineg_set:
            assert point_kept(point, ineg)

    print("Finished :" + "  {} inequalities".format(len(ineg_set)))

    # pure python results
    n = int(n)
    pure_ddt = dict()
    for a in range(int(1) << int(in_size)):
        for b in range(int(1) << int(out_size)):
            pure_ddt[int(a), int(b)] = int(ddt[a, b])
    in_size = int(in_size)
    out_size = int(out_size)

    return (
        in_size,
        out_size,
        pure_ddt,
        ineg_set,
        point_set,
        ineg_to_points,
        point_to_inegs,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute big set of inequalities with convex hull technique."
    )
    parser.add_argument(
        "sbox_name", type=str,
    )
    parser.add_argument(
        "nb_faces", type=int,
    )
    args = parser.parse_args()

    res = inequalities(sbox(args.sbox_name), args.nb_faces)

    print("Writing pickle file.")
    with open("{}_hull_{}.pkl".format(args.sbox_name, args.nb_faces), "wb") as f:
        pickle.dump(res, f, int(3))
