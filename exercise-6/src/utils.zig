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

pub fn InstanceCollection(comptime TypeCollection: type) type {
    return struct {
        pub fn call(_: @This(), k: Index, x: Scalar) Scalar {
            return TypeCollection.call(k, x);
        }
    };
}

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

pub fn OrtogonalCollection(comptime BaseCollection: type) type {
    return struct {
        const Self = @This();

        collection_: BaseCollection,
        order_: Index,
        coeficients_: []Scalar,
        diagonal_: []Scalar,

        pub fn init(allocator: std.mem.Allocator, order: Index) !Self {
            return Self{
                .collection_ = undefined,
                .order_ = order,
                .coeficients_ = try allocator.alloc(
                    Scalar,
                    order * (order - 1) / 2,
                ),
                .diagonal_ = try allocator.alloc(Scalar, order),
            };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.coeficients_);
            allocator.free(self.diagonal_);
        }

        pub fn compute(
            self: *Self,
            collection: BaseCollection,
            points: []Scalar,
        ) void {
            self.collection_ = collection;

            self.diagonal_[0] = 0;

            for (0..points.len) |j| {
                const v = self.collection_.call(0, points[j]);
                self.diagonal_[0] += v * v;
            }

            for (1..self.order_) |k| {
                const ks = k * (k - 1) / 2;

                for (0..k) |i| {
                    self.coeficients_[ks + i] = 0;
                }

                for (0..k) |i| {
                    var c: Scalar = 0;

                    for (0..points.len) |j| {
                        c += self.collection_.call(k, points[j]) *
                            self.call(i, points[j]);
                    }

                    self.coeficients_[ks + i] = -c / self.diagonal_[i];
                }

                for (0..k) |i| {
                    for (i + 1..k) |j| {
                        self.coeficients_[ks + i] +=
                            self.coeficients_[ks + j] *
                            self.coeficients_[(j * (j -% 1) / 2) + i];
                    }
                }

                self.diagonal_[k] = 0;

                for (0..points.len) |j| {
                    const v = self.call(k, points[j]);
                    self.diagonal_[k] += v * v;
                }
            }

            std.debug.print("{any}\n", .{self});
        }

        pub fn call(self: Self, k: Index, x: Scalar) Scalar {
            std.debug.assert(k < self.order_);
            var r: Scalar = self.collection_.call(k, x);
            const ks = k * (k -% 1) / 2;

            for (0..k) |i| {
                r += self.coeficients_[ks + i] * self.collection_.call(i, x);
            }

            return r;
        }

        pub fn diagonal(self: Self, k: Index) Scalar {
            return self.diagonal_[k];
        }
    };
}
