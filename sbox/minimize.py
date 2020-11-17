import argparse
import pickle
from gurobipy import *


def build_model(ineq_set, point_set, ineq_to_points, point_to_ineqs, start=None):
    """
    For the mode milp, this function builds the model as explained in ST17a.
    """
    model = Model("small milp")
    # model.hideOutput()

    z = {}
    for ineq in ineq_set:
        z[ineq] = model.addVar(name="z_{}".format(ineq), vtype=GRB.BINARY)
    if start is not None:
        for ineq in ineq_set:
            if ineq in start:
                z[ineq].start = 1.0
            else:
                z[ineq].start = 0.0

    for point in point_set:
        model.addConstr(
            quicksum(z[ineq] for ineq in point_to_ineqs[point]) >= 1,
            name="point {}".format(point),
        )
    return (model, z)


def optimize(
    ineq_set,
    point_set,
    ineq_to_points,
    point_to_ineqs,
    number=None,
    start=None,
    threshold=None,
):
    """
    Main function for the mode milp.
    ineq_set: big set of inequalities
    point_set: set of impossible points
    ineq_to_points: the map from inequalities to the impossible points they discard
    point_to_ineqs: the map from impossible points to the inequalities that discard them
    number:
        None -> minimize as much as possible
        _ -> the exact number of inequalities desired in the sbox model
    start: optional starting solution (usually given by the greedy algorithm)
    threshold: optional parameter to speed up the computation when ineq_set is too large
        only keeps the <threshold> best inequalities for each impossible point for
        the computation
    """

    # If there is a threshold, we update ineq_set and point_to_ineqs
    if threshold is not None:
        mut_ineq_set = set()
        mut_point_to_ineqs = {}

        # For each point, we keep the <threshold> best ineqs
        for point in point_set:
            key = lambda ineq: len(ineq_to_points[ineq])
            sorted_ineqs = sorted(list(point_to_ineqs[point]), key=key)
            best_ineqs = sorted_ineqs[:threshold]
            mut_point_to_ineqs[point] = set(best_ineqs)
            mut_ineq_set.update(set(best_ineqs))

        # And we keep the ineqs in the start set
        if start is not None:
            for ineq in start:
                mut_ineq_set.add(ineq)
                for point in ineq_to_points[ineq]:
                    mut_point_to_ineqs[point].add(ineq)

        ineq_set = mut_ineq_set
        point_to_ineqs = mut_point_to_ineqs

    model, z = build_model(ineq_set, point_set, ineq_to_points, point_to_ineqs, start)

    if number is None:
        model.setObjective(quicksum(z[ineq] for ineq in ineq_set), GRB.MINIMIZE)
    else:
        model.addConstr(quicksum(z[ineq] for ineq in ineq_set) <= number)
        model.setObjective(
            quicksum(z[ineq] * len(ineq_to_points[ineq]) for ineq in ineq_set),
            GRB.MAXIMIZE,
        )

    model.optimize()

    status = model.getAttr(GRB.Attr.Status)

    if status == GRB.INF_OR_UNBD or status == GRB.INFEASIBLE or status == GRB.UNBOUNDED:
        print("The model cannot be solved because it is infeasible or unbounded")
        sys.exit(1)
    if status != GRB.OPTIMAL:
        print("Optimization was stopped with status ", status)
        sys.exit(1)

    final_ineq_set = set()
    for ineq in ineq_set:
        if z[ineq].x >= 0.5:
            final_ineq_set.add(ineq)

    return final_ineq_set


def greedy_start(
    ineq_set, point_set, ineq_to_points,
):
    """
    Main function for the mode greedy.
    See optilize() for parameters description.
    """

    mut_point_set = point_set.copy()
    mut_ineq_set = ineq_set.copy()
    mut_ineq_to_points = ineq_to_points.copy()

    greedy_ineqs = set()

    while len(mut_point_set) > 0:
        # Search for best ineq
        best = mut_ineq_set.pop()
        mut_ineq_set.add(best)
        for ineq in mut_ineq_set:
            if len(mut_ineq_to_points[ineq]) > len(mut_ineq_to_points[best]):
                best = ineq
        points_of_best = mut_ineq_to_points.pop(best)
        mut_ineq_set.remove(best)
        greedy_ineqs.add(best)

        # Remove the points discarded by best from point sets
        mut_point_set.difference_update(points_of_best)
        for ineq in mut_ineq_set:
            mut_ineq_to_points[ineq].difference_update(points_of_best)

    return greedy_ineqs


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Reduction from big sets of inequalities to sbox models."
    )
    parser.add_argument(
        "ineq_file",
        type=str,
        help="Pickle file containing tuple "
        + "(in_size, out_size, ddt, ineq_set, point_set, "
        + "ineq_to_points, point_to_ineqs).",
    )
    parser.add_argument(
        "-m",
        type=str,
        dest="mode",
        default="greedy",
        choices=["milp", "greedy"],
        help="Either milp or greedy.",
    )
    parser.add_argument(
        "-n",
        type=int,
        dest="number",
        help="Number of inequalities if the chosen mode is milp.",
    )
    parser.add_argument(
        "-t",
        type=int,
        dest="threshold",
        help="Number of inequalities per point kept for building the model"
        + " if the chosen mode is milp.",
    )
    parser.add_argument(
        "-sf",
        type=str,
        dest="start_file",
        help="Pickle file with a starting solution (set of inequalities)"
        + " if the chosen mode is milp.",
    )
    args = parser.parse_args()

    output_file = "sbox_{}".format(args.ineq_file)

    with open(args.ineq_file, "rb") as f:
        (
            in_size,
            out_size,
            ddt,
            ineq_set,
            point_set,
            ineq_to_points,
            point_to_ineqs,
        ) = pickle.load(f)

    if args.mode == "greedy":
        final_set = greedy_start(ineq_set, point_set, ineq_to_points)
        output_file = "greedy_" + output_file
    else:
        if args.start_file is not None:
            with open(args.start_file, "rb") as f:
                (_, _, ddt_start, ineq_start) = pickle.load(f)
            assert ddt_start == ddt
        else:
            ineq_start = None

        final_set = optimize(
            ineq_set,
            point_set,
            ineq_to_points,
            point_to_ineqs,
            number=args.number,
            start=ineq_start,
            threshold=args.threshold,
        )

        if args.number is None:
            output_file = "milp_minimal_" + output_file
        else:
            output_file = "milp_{}_".format(args.number) + output_file

    with open(output_file, "wb") as f:
        pickle.dump((in_size, out_size, ddt, final_set), f, 3)
