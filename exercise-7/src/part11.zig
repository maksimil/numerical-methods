const std = @import("std");
const config = @import("config.zig");
const utils = @import("utils.zig");

const Scalar = config.Scalar;

pub fn Run() !void {
    var output_dir = try std.fs.cwd().makeOpenPath("output", .{});
    const file = try output_dir.createFile("part11.csv", .{});
    defer file.close();
    const writer = file.writer();

    try config.stdout.print("\x1B[1;32m>>> Part 1.1\x1B[0m\n", .{});

    try writer.print("Evaluations,LTrianglesErr,RTrianglesErr," ++
        "MTrianglesErr,TrapezoidErr,SimpsonErr\n", .{});

    const h = config.kTaskB - config.kTaskA;

    for (0..config.kSections.len) |ni| {
        const n = config.kSections[ni];

        const left_triangles = utils.ChunkedSum(struct {
            n: usize,
            pub fn call(self: @This(), k: usize) Scalar {
                return config.TaskF(
                    config.kTaskA +
                        @as(Scalar, @floatFromInt(k)) /
                        @as(Scalar, @floatFromInt(self.n)) * h,
                );
            }
        }, .{ .n = n }, 0, n, utils.kChunkSize) *
            h / @as(Scalar, @floatFromInt(n));

        const right_triangles = utils.ChunkedSum(struct {
            n: usize,
            pub fn call(self: @This(), k: usize) Scalar {
                return config.TaskF(
                    config.kTaskA +
                        @as(Scalar, @floatFromInt(k + 1)) /
                        @as(Scalar, @floatFromInt(self.n)) * h,
                );
            }
        }, .{ .n = n }, 0, n, utils.kChunkSize) *
            h / @as(Scalar, @floatFromInt(n));

        const middle_triangles = utils.ChunkedSum(struct {
            n: usize,
            pub fn call(self: @This(), k: usize) Scalar {
                return config.TaskF(
                    config.kTaskA +
                        (@as(Scalar, @floatFromInt(k)) + 0.5) /
                        @as(Scalar, @floatFromInt(self.n)) * h,
                );
            }
        }, .{ .n = n }, 0, n, utils.kChunkSize) *
            h / @as(Scalar, @floatFromInt(n));

        const trapezoid = (2.0 * utils.ChunkedSum(struct {
            n: usize,
            pub fn call(self: @This(), k: usize) Scalar {
                return config.TaskF(
                    config.kTaskA +
                        @as(Scalar, @floatFromInt(k)) /
                        @as(Scalar, @floatFromInt(self.n)) * h,
                );
            }
        }, .{ .n = n }, 1, n - 1, utils.kChunkSize) +
            config.TaskF(config.kTaskA) + config.TaskF(config.kTaskB)) /
            2.0 * h / @as(Scalar, @floatFromInt(n));

        const simpson = (trapezoid + 2.0 * middle_triangles) / 3.0;

        try writer.print("{d},{e},{e},{e},,\n", .{
            n,
            @abs(config.kTaskAnswer - left_triangles),
            @abs(config.kTaskAnswer - right_triangles),
            @abs(config.kTaskAnswer - middle_triangles),
        });

        try writer.print("{d},,,,{e},\n", .{
            n + 1,
            @abs(config.kTaskAnswer - trapezoid),
        });

        try writer.print("{d},,,,,{e}\n", .{
            2 * n + 1,
            @abs(config.kTaskAnswer - simpson),
        });

        try config.stdout.print(
            "Computed N={d:10}, Ni={d:3}/{d}\n",
            .{
                n,
                ni + 1,
                config.kSections.len,
            },
        );
    }

    try config.stdout.print("\x1B[34mWrote to output/part11.csv\x1B[0m\n", .{});
    try config.stdout.print("\x1B[1;32m<<< End\x1B[0m\n\n", .{});
}
