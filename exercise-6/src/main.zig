const std = @import("std");
const config = @import("config.zig");
const matrix = @import("matrix.zig");
const lu = @import("lu.zig");

const Scalar = f64;

pub fn main() !void {
    var arena_allocator = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena_allocator.deinit();

    const allocator = arena_allocator.allocator();

    var m = try matrix.RowMajorMatrix.init(allocator, 2);
    defer m.deinit(allocator);

    m.set_zero();
    m.at_mut(0, 0).* = 1;
    m.at_mut(1, 0).* = 2;
    m.at_mut(1, 1).* = 1;

    std.debug.print("{any}\n", .{m});

    var d = try lu.LUDecomposition.init(allocator, 2);
    _ = d.factorize(matrix.RowMajorMatrix, m);

    std.debug.print("{any}\n", .{d});

    defer d.deinit(allocator);
}
