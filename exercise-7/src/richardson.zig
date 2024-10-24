const std = @import("std");
const config = @import("config.zig");
const lu = @import("lu/lu.zig");
const mt = @import("lu/matrix.zig");

const Scalar = config.Scalar;

fn powi(x: Scalar, n: usize) Scalar {
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

pub fn RichardsonExtrapolation(
    steps: []const Scalar,
    values: []const Scalar,
    convergence: usize,
) !Scalar {
    const r = steps.len;
    std.debug.assert(values.len == r);

    if (r == 1) {
        return config.kNan;
    }

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
        try config.stdout.print("Lowered extrapolation r={d:3}\n", .{r - 1});
        return RichardsonExtrapolation(steps[1..], values[1..], convergence);
    }

    decomposition.solve(mt.Vector, &vector);

    return vector.at(r - 1);
}

const kMulOrder = 2;
const kMulOrderScalar = @as(Scalar, @floatFromInt(kMulOrder));
const kLogMulOrder = @log(kMulOrderScalar);

pub fn ProcessFormula(
    comptime Formula: type,
    target_eps: Scalar,
) !void {
    var steps = std.ArrayList(Scalar).init(config.allocator);
    defer steps.deinit();

    var values = std.ArrayList(Scalar).init(config.allocator);
    defer values.deinit();

    const h = config.kTaskB - config.kTaskA;

    try steps.append(h);

    try values.append(Formula.call(1));

    var n: usize = 1;
    var err: Scalar = 2 * target_eps;

    while (@abs(err) >= target_eps) {
        n *= kMulOrder;

        try steps.append(h / @as(Scalar, @floatFromInt(n)));
        try values.append(Formula.call(n));

        const extrapolation_value =
            try RichardsonExtrapolation(
            steps.items,
            values.items,
            Formula.kConvergence,
        );
        err = extrapolation_value - values.items[values.items.len - 1];

        const eitkin: Scalar = blk: {
            if (steps.items.len < 3) break :blk config.kNan;

            const s1 = values.items[values.items.len - 3];
            const s2 = values.items[values.items.len - 2];
            const s3 = values.items[values.items.len - 1];

            const v = -@log(@abs((s3 - s2) / (s2 - s1))) / kLogMulOrder;

            if ((s3 - s2) / (s2 - s1) < 0) {
                break :blk -v;
            } else {
                break :blk v;
            }
        };

        try config.stdout.print(
            "r={d:3} n={d:6} err={e:10.3} true_err={e:10.3} " ++
                "ep_err={e:10.3} eitkin={e:10.3}\n",
            .{
                steps.items.len,
                n,
                err,
                config.kTaskAnswerP - values.items[values.items.len - 1],
                config.kTaskAnswerP - extrapolation_value,
                eitkin,
            },
        );
    }

    try config.stdout.print(
        "Computed, h={e:12.5}\n",
        .{steps.items[steps.items.len - 1]},
    );
}
