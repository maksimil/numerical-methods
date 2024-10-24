const config = @import("config.zig");
const std = @import("std");
const lu = @import("lu/lu.zig");

const Scalar = config.Scalar;

pub const kChunkSize = 32;
// pub const kChunkSize = std.math.maxInt(usize);

pub fn ChunkedSum(
    comptime GeneratorType: type,
    generator: GeneratorType,
    start: usize,
    n: usize,
    chunksize: usize,
) Scalar {
    var sum: Scalar = 0.0;
    if (n <= chunksize) {
        for (0..n) |i| {
            sum += generator.call(start + i);
        }
    } else {
        const step = n / chunksize;
        for (0..(chunksize - 1)) |chunk| {
            sum += ChunkedSum(
                GeneratorType,
                generator,
                start + step * chunk,
                step,
                chunksize,
            );
        }
        sum += ChunkedSum(
            GeneratorType,
            generator,
            start + step * (chunksize - 1),
            n - step * (chunksize - 1),
            chunksize,
        );
    }
    return sum;
}

// Solve
//     [a0 a1 a2; a3 a4 a5; a6 a7 a8] [x0; x1; x2] = [b0; b1; b2]
// where x is both b and then x

pub fn Solve3(a: []const Scalar, x: []Scalar) bool {
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

// Indefinite integrals of x^k p(x)
// all the functions assume alpha=0

fn root6(x: Scalar) Scalar {
    return std.math.cbrt(@sqrt(x));
}

const b = config.kTaskB;

pub fn IndInt0(x: Scalar) Scalar {
    return -6.0 * root6(b - x);
}

pub fn IndInt1(x: Scalar) Scalar {
    return -6.0 / 7.0 * root6(b - x) * (6.0 * b + x);
}

pub fn IndInt2(x: Scalar) Scalar {
    return -6.0 / 91.0 * root6(b - x) *
        (72.0 * b * b + 12.0 * b * x + 7.0 * x * x);
}

pub fn IndInt3(x: Scalar) Scalar {
    return -6.0 / 1729.0 * root6(b - x) *
        (1296.0 * b * b * b + 216.0 * b * b * x +
        126.0 * b * x * x + 91.0 * x * x * x);
}

pub fn IndInt4(x: Scalar) Scalar {
    return -6.0 / 43225.0 * root6(b - x) *
        (31104.0 * b * b * b * b + 5184.0 * b * b * b * x +
        3024.0 * b * b * x * x + 2184.0 * b * x * x * x +
        1729.0 * x * x * x * x);
}

pub fn IndInt5(x: Scalar) Scalar {
    const p1 = root6(b - x);
    const p2 = p1 * p1;
    const p4 = p2 * p2;
    const p8 = p4 * p4;
    const p16 = p8 * p8;
    return (-6.0 * (b * b * b * b * b) * p1) +
        (30.0 / 7.0 * (b * b * b * b) * p4 * p2 * p1) +
        (-60.0 / 13.0 * (b * b * b) * p8 * p4 * p1) +
        (60.0 / 19.0 * (b * b) * p16 * p2 * p1) +
        (-6.0 / 5.0 * b * p16 * p8 * p1) +
        (6.0 / 31.0 * p16 * p8 * p4 * p2 * p1);
}
