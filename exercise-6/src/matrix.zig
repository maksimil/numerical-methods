const std = @import("std");
const config = @import("config.zig");

const Scalar = config.Scalar;
const Index = config.Index;

pub const Vector = struct {
    const Self = @This();

    data_: []Scalar,

    pub fn from_slice(data: []Scalar) Self {
        return Self{ .data_ = data };
    }

    pub fn init(allocator: std.mem.Allocator, dim: Index) !Self {
        return from_slice(try allocator.alloc(Scalar, dim));
    }

    pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        allocator.free(self.data_);
    }

    pub fn dimension(self: Self) Index {
        return self.data_.len;
    }

    pub fn set_zero(self: *Self) void {
        @memset(self.data_[0..self.dimension()], 0);
    }

    pub fn at(self: Self, i: Index) Scalar {
        return self.data_[i];
    }

    pub fn at_mut(self: *Self, i: Index) *Scalar {
        return &self.data_[i];
    }
};

pub const RowMajorMatrix = struct {
    const Self = @This();

    data_: []Scalar,
    dimension_: Index,

    pub fn init(allocator: std.mem.Allocator, dim: Index) !Self {
        const data = try allocator.alloc(Scalar, dim * dim);
        return Self{ .data_ = data, .dimension_ = dim };
    }

    pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        allocator.free(self.data_);
    }

    pub fn set_zero(self: *Self) void {
        @memset(self.data_[0..(self.dimension() * self.dimension())], 0);
    }

    pub fn dimension(self: Self) Index {
        return self.dimension_;
    }

    pub fn at(self: Self, row: Index, col: Index) Scalar {
        return self.data_[row * self.dimension() + col];
    }

    pub fn at_mut(self: *Self, row: Index, col: Index) *Scalar {
        return &self.data_[row * self.dimension() + col];
    }
};

pub fn RowOf(comptime Storage: type) type {
    return struct {
        const Self = @This();

        matrix: Storage,
        row: Index,

        pub fn init(matrix: Storage, row: Index) Self {
            return Self{ .matrix = matrix, .row = row };
        }

        pub fn dimension(self: Self) Index {
            return self.matrix.dimension();
        }

        pub fn at(self: Self, i: Index) Scalar {
            return self.matrix.at(self.row, i);
        }

        pub fn at_mut(self: *Self, i: Index) *Scalar {
            return self.matrix.at_mut(self.row, i);
        }
    };
}

pub fn ColOf(comptime Storage: type) type {
    return struct {
        const Self = @This();

        matrix: Storage,
        col: Index,

        pub fn init(matrix: Storage, col: Index) Self {
            return Self{ .matrix = matrix, .col = col };
        }

        pub fn dimension(self: Self) Index {
            return self.matrix.dimension();
        }

        pub fn at(self: Self, i: Index) Scalar {
            return self.matrix.at(i, self.col);
        }

        pub fn at_mut(self: *Self, i: Index) *Scalar {
            return self.matrix.at_mut(i, self.col);
        }
    };
}

pub const Permutation = struct {
    const Self = @This();

    forward_: []Index,
    inverse_: []Index,
    dimension_: Index,

    pub fn init(allocator: std.mem.Allocator, dimension: Index) !Self {
        return Self{
            .forward_ = try allocator.alloc(Index, dimension),
            .inverse_ = try allocator.alloc(Index, dimension),
            .dimension_ = dimension,
        };
    }

    pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        allocator.free(self.forward_);
        allocator.free(self.inverse_);
    }

    pub fn set_identity(self: *Self) void {
        for (0..self.dimension_) |i| {
            self.forward_[i] = i;
            self.inverse_[i] = i;
        }
    }

    pub fn set(self: *Self, i: Index, pi: Index) void {
        self.forward_[i] = pi;
        self.inverse_[pi] = i;
    }

    pub fn consistent_set(self: *Self, i: Index, pi: Index) void {
        self.swap(i, self.inv(pi));
    }

    pub fn swap(self: *Self, i: Index, j: Index) void {
        const pi = self.perm(i);
        const pj = self.perm(j);
        self.set(i, pj);
        self.set(j, pi);
    }

    pub fn perm(self: Self, i: Index) Index {
        return self.forward_[i];
    }

    pub fn inv(self: Self, pi: Index) Index {
        return self.inverse_[pi];
    }
};
