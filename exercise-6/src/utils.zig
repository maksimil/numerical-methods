const config = @import("config.zig");
const matrix = @import("matrix.zig");
const std = @import("std");

const Scalar = config.Scalar;
const Index = config.Index;

pub const Triangularity = enum { Upper, Lower };

pub fn dot_product(
    comptime A: type,
    comptime B: type,
    a: A,
    b: B,
) Scalar {
    std.debug.assert(a.dimension() == b.dimension());

    var r: Scalar = 0;

    for (0..a.dimension()) |i| {
        r += a.at(i) * b.at(i);
    }

    return r;
}

pub fn triangular_solve(
    comptime MatrixRefType: type,
    comptime VectorRefType: type,
    comptime triangularity: Triangularity,
    mtx: MatrixRefType,
    permutation: matrix.Permutation,
    vector: *VectorRefType,
) void {
    const dimension = mtx.dimension();

    for (0..dimension) |kk| {
        var k: Index = undefined;

        if (triangularity == Triangularity.Upper) {
            k = permutation.perm(dimension - 1 - kk);
        } else {
            k = permutation.perm(kk);
        }

        const v = vector.at(k);
        vector.at_mut(k).* = 0;
        const dot = dot_product(
            matrix.RowOf(matrix.RowMajorMatrix),
            VectorRefType,
            matrix.RowOf(matrix.RowMajorMatrix).init(mtx, k),
            vector.*,
        );
        vector.at_mut(k).* = (v - dot) / mtx.at(k, k);
    }
}

pub fn powi(x: Scalar, n: Index) Scalar {
    var x2k: Scalar = x; // x^{2k};
    var exp: Index = n;

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

pub fn FnType(comptime f: fn (Scalar) Scalar) type {
    return struct {
        pub fn call(x: Scalar) Scalar {
            return f(x);
        }
    };
}

pub const PolynomialCollection = struct {
    pub fn call(k: Index, x: Scalar) Scalar {
        return powi(x, k);
    }
};

pub fn rand_scalar(rnd: std.Random, a: Scalar, b: Scalar) Scalar {
    return a + (b - a) * rnd.float(Scalar);
}

pub fn collection_call(
    comptime CollectionType: type,
    collection: CollectionType,
    coefs: []const Scalar,
    x: Scalar,
) Scalar {
    var r: Scalar = 0;

    for (0..coefs.len) |k| {
        r += coefs[k] * collection.call(k, x);
    }

    return r;
}
