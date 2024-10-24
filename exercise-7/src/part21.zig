const std = @import("std");
const config = @import("config.zig");
const lu = @import("lu/lu.zig");
const mt = @import("lu/matrix.zig");
const weighted = @import("weighted.zig");

const Scalar = config.Scalar;

const TARGET_EPS = 1e-6;
const NC_CONVERGENCE = 3;
const GAUSS_CONVERGENCE = 5;

pub fn powi(x: Scalar, n: usize) Scalar {
    var x2k = x; // x^{2k};
    var exp = n;

    var r: Scalar = 1;

    while (exp > 0) {
        if (exp & 1 == 1) {
            r *= x2k;
        }

        x2k *= x2k;

        exp >>= 1;
    }

    return r;
}

fn RichardsonExtrapolation(
    steps: []const Scalar,
    values: []const Scalar,
    convergence: usize,
) !Scalar {
    const r = steps.len;
    std.debug.assert(values.len == r);

    const matrix_data = try config.allocator.alloc(Scalar, r * r);
    defer config.allocator.free(matrix_data);

    const vector_data = try config.allocator.alloc(Scalar, r);
    defer config.allocator.free(vector_data);

    var matrix = mt.RowMajorMatrix.from_slice(matrix_data, r);
    var vector = mt.Vector.from_slice(vector_data);

    for (0..r) |i| {
        for (0..r - 1) |j| {
            matrix.at_mut(i, j).* = powi(steps[i], convergence + j);
        }

        matrix.at_mut(i, r - 1).* = -1.0;

        vector.at_mut(i).* = -values[i];
    }

    var decomposition = try lu.LUDecomposition.init(config.allocator, r);
    defer decomposition.deinit(config.allocator);

    const factored = decomposition.factorize(mt.RowMajorMatrix, matrix);

    if (!factored) {
        return config.kNan;
    }

    decomposition.solve(mt.Vector, &vector);

    return vector.at(r - 1);
}

pub fn Run() !void {
    try config.stdout.print("\x1B[1;32m>>> Part 2.1\x1B[0m\n", .{});

    var steps = std.ArrayList(Scalar).init(config.allocator);
    defer steps.deinit();

    var values = std.ArrayList(Scalar).init(config.allocator);
    defer values.deinit();

    const h = config.kTaskB - config.kTaskA;

    try steps.append(h);
    try steps.append(h / 2.0);

    try values.append(weighted.ComputeNC(1));
    try values.append(weighted.ComputeNC(2));

    var n: usize = 2;
    var err = try RichardsonExtrapolation(
        steps.items,
        values.items,
        NC_CONVERGENCE,
    ) - values.items[1];

    while (@abs(err) >= TARGET_EPS) {
        try config.stdout.print(
            "r={d:3} n={d:6} err={e:10.3} true_err={e:10.3}\n",
            .{
                steps.items.len,
                n,
                err,
                config.kTaskAnswerP - values.items[values.items.len - 1],
            },
        );

        n *= 2;

        try steps.append(h / @as(Scalar, @floatFromInt(n)));
        try values.append(weighted.ComputeNC(n));

        err = try RichardsonExtrapolation(steps.items, values.items, NC_CONVERGENCE) -
            values.items[values.items.len - 1];
    }

    try config.stdout.print(
        "r={d:3} n={d:6} err={e:10.3} true_err={e:10.3}\n",
        .{
            steps.items.len,
            n,
            err,
            config.kTaskAnswerP - values.items[values.items.len - 1],
        },
    );

    try config.stdout.print(
        "Computed, h={e:12.5}\n",
        .{steps.items[steps.items.len - 1]},
    );

    try config.stdout.print("\x1B[1;32m<<< End\x1B[0m\n\n", .{});
}
