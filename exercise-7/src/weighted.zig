const std = @import("std");
const config = @import("config.zig");
const utils = @import("utils.zig");

const Scalar = config.Scalar;

// Solve
//     [a0 a1 a2; a3 a4 a5; a6 a7 a8] [x0; x1; x2] = [b0; b1; b2]
// where x is both b and then x
fn Solve3(a: []const Scalar, x: []Scalar) bool {
    var matrix: [9]Scalar = undefined;
    var p0: usize = undefined;
    var p1: usize = undefined;
    var p2: usize = undefined;
    var coef: Scalar = undefined;

    @memcpy(&matrix, a);

    // Gaussian elimination
    {
        p0 = 0;

        for (1..3) |k| {
            if (@abs(matrix[0 + 3 * p0]) < @abs(matrix[0 + 3 * k])) {
                p0 = k;
            }
        }

        if (@abs(matrix[0 + 3 * p0]) < config.MIN_PIVOT) {
            return false;
        }

        coef = undefined;

        for (0..3) |k| {
            if (p0 != k) {
                coef = matrix[0 + 3 * k] / matrix[0 + 3 * p0];

                // matrix[0 + 3 * k] = 0;
                matrix[1 + 3 * k] -= coef * matrix[1 + 3 * p0];
                matrix[2 + 3 * k] -= coef * matrix[2 + 3 * p0];
                x[k] -= coef * x[p0];
            }
        }
    }

    {
        if (@abs(matrix[1 + 3 * ((p0 + 1) % 3)]) <
            @abs(matrix[1 + 3 * ((p0 + 2) % 3)]))
        {
            p1 = (p0 + 2) % 3;
        } else {
            p1 = (p0 + 1) % 3;
        }

        if (@abs(matrix[1 + 3 * p1]) < config.MIN_PIVOT) {
            return false;
        }

        p2 = 3 - p1 - p0;

        coef = matrix[1 + 3 * p2] / matrix[1 + 3 * p1];

        // matrix[1 + 3 * p2] = 0;
        matrix[2 + 3 * p2] -= coef * matrix[2 + 3 * p1];
        x[p2] -= coef * x[p1];
    }

    {
        if (@abs(matrix[2 + 3 * p2]) < config.MIN_PIVOT) {
            return false;
        }

        std.debug.assert(@abs(matrix[0 + 3 * p0]) >= config.MIN_PIVOT);
        std.debug.assert(@abs(matrix[1 + 3 * p1]) >= config.MIN_PIVOT);
        std.debug.assert(@abs(matrix[2 + 3 * p2]) >= config.MIN_PIVOT);

        const r: [3]Scalar = .{ x[0], x[1], x[2] };
        x[2] = (r[p2]) / matrix[2 + 3 * p2];
        x[1] = (r[p1] - matrix[2 + 3 * p1] * x[2]) / matrix[1 + 3 * p1];
        x[0] = (r[p0] - matrix[1 + 3 * p0] * x[1] -
            matrix[2 + 3 * p0] * x[2]) / matrix[0 + 3 * p0];
    }

    return true;
}

fn CardanoFormula(a0: Scalar, a1: Scalar, a2: Scalar, roots: []Scalar) void {
    const p = (3.0 * a1 - a2 * a2) / 9.0;
    const q = (a2 * a2 * a2) / 27.0 - a2 * a1 / 6.0 + a0 / 2.0;

    var r: Scalar = @sqrt(@abs(p));

    if (q < 0.0) {
        r = -r;
    }

    var v = q / (r * r * r);

    if (v >= 1.0) {
        v = 1.0;
    } else if (v <= -1.0) {
        v = -1.0;
    }

    const phi = std.math.acos(v);

    std.debug.assert(phi >= 0);
    std.debug.assert(phi <= std.math.pi);

    roots[0] = -a2 / 3.0 + 2.0 * r * @cos(phi / 3.0 + std.math.pi);
    roots[1] = -a2 / 3.0 + 2.0 * r * @cos(phi / 3.0 + std.math.pi / 3.0);
    roots[2] = -a2 / 3.0 + 2.0 * r * @cos(phi / 3.0 - std.math.pi / 3.0);

    if (r < 0) {
        std.mem.swap(Scalar, &roots[0], &roots[2]);
    }

    std.debug.assert(roots[0] <= roots[1]);
    std.debug.assert(roots[1] <= roots[2]);

    // config.stdout.print("{any}\n", .{roots}) catch unreachable;
}

fn QFCoefs(
    x: []const Scalar,
    coefs: []Scalar,
    start: Scalar,
    end: Scalar,
) bool {
    const matrix: [9]Scalar = .{
        1.0,
        1.0,
        1.0,
        x[0],
        x[1],
        x[2],
        x[0] * x[0],
        x[1] * x[1],
        x[2] * x[2],
    };

    coefs[0] = utils.IndInt0(end) - utils.IndInt0(start);
    coefs[1] = utils.IndInt1(end) - utils.IndInt1(start);
    coefs[2] = utils.IndInt2(end) - utils.IndInt2(start);

    const r = utils.Solve3(&matrix, coefs);

    return r;
}

fn ComputeByPoints(
    x0: Scalar,
    x1: Scalar,
    x2: Scalar,
    z0: Scalar,
    z1: Scalar,
) Scalar {
    var coefs: [3]Scalar = undefined;

    const matrix: [9]Scalar = .{
        1.0,     1.0,     1.0,
        x0,      x1,      x2,
        x0 * x0, x1 * x1, x2 * x2,
    };

    coefs[0] = utils.IndInt0(z1) - utils.IndInt0(z0);
    coefs[1] = utils.IndInt1(z1) - utils.IndInt1(z0);
    coefs[2] = utils.IndInt2(z1) - utils.IndInt2(z0);

    const r = Solve3(&matrix, &coefs);

    if (!r) return config.kNan;

    return coefs[0] * config.TaskF(x0) +
        coefs[1] * config.TaskF(x1) +
        coefs[2] * config.TaskF(x2);
}

pub fn NewtonCotes(z0: Scalar, z1: Scalar) Scalar {
    const x1 = (z0 + z1) / 2.0;
    return ComputeByPoints(z0, x1, z1, z0, z1);
}

pub fn Gauss(z0: Scalar, z1: Scalar) Scalar {
    const mu0 = utils.IndInt0(z1) - utils.IndInt0(z0);
    const mu1 = utils.IndInt1(z1) - utils.IndInt1(z0);
    const mu2 = utils.IndInt2(z1) - utils.IndInt2(z0);
    const mu3 = utils.IndInt3(z1) - utils.IndInt3(z0);
    const mu4 = utils.IndInt4(z1) - utils.IndInt4(z0);
    const mu5 = utils.IndInt5(z1) - utils.IndInt5(z0);

    const matrix: [9]Scalar = .{
        mu0, mu1, mu2,
        mu1, mu2, mu3,
        mu2, mu3, mu4,
    };

    var coefs: [3]Scalar = .{ -mu3, -mu4, -mu5 };

    const r = Solve3(&matrix, &coefs);

    if (!r) return config.kNan;

    var x: [3]Scalar = undefined;
    CardanoFormula(coefs[0], coefs[1], coefs[2], &x);

    return ComputeByPoints(x[0], x[1], x[2], z0, z1);
}

pub fn ComputeNC(n: usize) Scalar {
    const h = config.kTaskB - config.kTaskA;

    const generator = struct {
        n: usize,

        pub fn call(self: @This(), k: usize) Scalar {
            const z0 = config.kTaskA +
                @as(Scalar, @floatFromInt(k)) /
                @as(Scalar, @floatFromInt(self.n)) * h;

            const z1 = config.kTaskA +
                @as(Scalar, @floatFromInt(k + 1)) /
                @as(Scalar, @floatFromInt(self.n)) * h;

            return NewtonCotes(z0, z1);
        }
    };

    return utils.ChunkedSum(generator, .{ .n = n }, 0, n, utils.kChunkSize);
}

pub fn ComputeG(n: usize) Scalar {
    const h = config.kTaskB - config.kTaskA;

    const generator = struct {
        n: usize,

        pub fn call(self: @This(), k: usize) Scalar {
            const z0 = config.kTaskA +
                @as(Scalar, @floatFromInt(k)) /
                @as(Scalar, @floatFromInt(self.n)) * h;

            const z1 = config.kTaskA +
                @as(Scalar, @floatFromInt(k + 1)) /
                @as(Scalar, @floatFromInt(self.n)) * h;

            return Gauss(z0, z1);
        }
    };

    return utils.ChunkedSum(generator, .{ .n = n }, 0, n, utils.kChunkSize);
}
