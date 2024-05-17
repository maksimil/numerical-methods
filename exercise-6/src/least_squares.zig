const std = @import("std");
const config = @import("config.zig");
const matrix = @import("matrix.zig");
const lu = @import("lu.zig");

const Scalar = config.Scalar;

pub fn least_squares(
    comptime FunctionCollection: type,
    comptime Orthogonal: bool,
    allocator: std.mem.Allocator,
    coefs: []Scalar,
    functions: FunctionCollection,
    points: []const Scalar,
    data: []const Scalar,
) !bool {
    const order = coefs.len;

    // data vector fill
    for (0..order) |row| {
        coefs[row] = 0;

        for (0..points.len) |i| {
            coefs[row] +=
                functions.call(row, points[i]) * data[i];
        }
    }

    if (Orthogonal) {
        for (0..order) |k| {
            coefs[k] /= functions.diagonal(k);
        }
    } else {
        // gram matrix
        var gram_matrix = try matrix.RowMajorMatrix.init(allocator, order);
        defer gram_matrix.deinit(allocator);

        for (0..order) |row| {
            for (0..order) |col| {
                gram_matrix.at_mut(row, col).* = 0;

                for (0..points.len) |i| {
                    gram_matrix.at_mut(row, col).* +=
                        functions.call(row, points[i]) *
                        functions.call(col, points[i]);
                }
            }
        }

        // computation
        var decomposition = try lu.LUDecomposition.init(allocator, order);
        defer decomposition.deinit(allocator);

        const factored =
            decomposition.factorize(matrix.RowMajorMatrix, gram_matrix);

        if (!factored) {
            return false;
        }

        var coefs_slice = matrix.Vector.from_slice(coefs);
        decomposition.solve(matrix.Vector, &coefs_slice);
    }

    return true;
}
