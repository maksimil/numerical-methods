const config = @import("../config.zig");
const matrix = @import("matrix.zig");
const std = @import("std");

const Scalar = config.Scalar;

const Triangularity = enum { Upper, Lower };

fn dot_product(
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

fn triangular_solve(
    comptime MatrixRefType: type,
    comptime VectorRefType: type,
    comptime triangularity: Triangularity,
    mtx: MatrixRefType,
    permutation: matrix.Permutation,
    vector: *VectorRefType,
) void {
    const dimension = mtx.dimension();

    for (0..dimension) |kk| {
        var k: usize = undefined;

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

pub const LUDecomposition = struct {
    const Self = @This();

    const MatrixStorage = matrix.RowMajorMatrix;
    const WorkVector = matrix.Vector;

    lower_: MatrixStorage,
    upper_: MatrixStorage,
    permutation_: matrix.Permutation,

    work_row_used_: []bool,
    work_vector_: WorkVector,

    pub fn init(allocator: std.mem.Allocator, dimension: usize) !Self {
        return Self{
            .lower_ = try MatrixStorage.init(allocator, dimension),
            .upper_ = try MatrixStorage.init(allocator, dimension),
            .permutation_ = try matrix.Permutation.init(allocator, dimension),
            .work_row_used_ = try allocator.alloc(bool, dimension),
            .work_vector_ = try matrix.Vector.init(allocator, dimension),
        };
    }

    pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        self.lower_.deinit(allocator);
        self.upper_.deinit(allocator);
        self.permutation_.deinit(allocator);

        allocator.free(self.work_row_used_);
        self.work_vector_.deinit(allocator);
    }

    // true if computes, false otherwise
    pub fn factorize(
        self: *Self,
        comptime Matrix: type,
        base: Matrix,
    ) bool {
        const dimension = self.lower_.dimension();

        // memory stuff
        var row_used = self.work_row_used_;
        @memset(row_used[0..dimension], false);

        var work_vector = self.work_vector_;
        work_vector.set_zero();

        self.permutation_.set_identity();
        self.lower_.set_zero();
        self.upper_.set_zero();

        // setting L to identity
        for (0..dimension) |i| {
            self.lower_.at_mut(i, i).* = 1;
        }

        // computing LU
        for (0..dimension) |col| {
            for (0..dimension) |i| {
                work_vector.at_mut(i).* = base.at(i, col);
            }

            triangular_solve(
                MatrixStorage,
                matrix.Vector,
                Triangularity.Lower,
                self.lower_,
                self.permutation_,
                &work_vector,
            );

            var pivot_row: usize = 0;

            for (0..dimension) |i| {
                if (!row_used[i] and
                    (row_used[pivot_row] or
                    @abs(work_vector.at(i)) > @abs(work_vector.at(pivot_row))))
                {
                    pivot_row = i;
                }
            }

            const pivot_value = work_vector.at(pivot_row);

            if (@abs(pivot_value) < config.MIN_PIVOT) {
                return false;
            }

            row_used[pivot_row] = true;

            for (0..dimension) |i| {
                if (row_used[i]) {
                    self.upper_.at_mut(i, pivot_row).* =
                        work_vector.at(i);
                } else {
                    self.lower_.at_mut(i, pivot_row).* =
                        work_vector.at(i) / pivot_value;
                }
            }

            self.permutation_.consistent_set(col, pivot_row);
        }

        return true;
    }

    pub fn solve(
        self: *Self,
        comptime VectorType: type,
        vector: *VectorType,
    ) void {
        triangular_solve(
            MatrixStorage,
            matrix.Vector,
            Triangularity.Lower,
            self.lower_,
            self.permutation_,
            vector,
        );

        triangular_solve(
            MatrixStorage,
            matrix.Vector,
            Triangularity.Upper,
            self.upper_,
            self.permutation_,
            vector,
        );

        const dimension = self.lower_.dimension();

        for (0..dimension) |i| {
            self.work_vector_.at_mut(i).* = vector.at(self.permutation_.perm(i));
        }

        for (0..dimension) |i| {
            vector.at_mut(i).* = self.work_vector_.at(i);
        }
    }
};
