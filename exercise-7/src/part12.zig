const std = @import("std");
const config = @import("config.zig");
const utils = @import("utils.zig");
const weighted = @import("weighted.zig");

const Scalar = config.Scalar;

pub fn Run() !void {
    var output_dir = try std.fs.cwd().makeOpenPath("output", .{});
    const file = try output_dir.createFile("part12.csv", .{});
    defer file.close();
    const writer = file.writer();

    try config.stdout.print("\x1B[1;32m>>> Part 1.2\x1B[0m\n", .{});

    try writer.print("Evaluations,NewtonCotesErr,GaussErr,HybridErr\n", .{});

    const h = config.kTaskB - config.kTaskA;

    const hybrid_generator = struct {
        n: usize,

        pub fn call(self: @This(), k: usize) Scalar {
            const z0 = config.kTaskA +
                @as(Scalar, @floatFromInt(k)) /
                @as(Scalar, @floatFromInt(self.n)) * h;

            const z1 = config.kTaskA +
                @as(Scalar, @floatFromInt(k + 1)) /
                @as(Scalar, @floatFromInt(self.n)) * h;

            const g = weighted.Gauss(z0, z1);

            if (std.math.isNan(g)) {
                return weighted.NewtonCotes(z0, z1);
            } else {
                return g;
            }
        }
    };

    for (0..config.kSections.len) |ni| {
        const n = config.kSections[ni];

        const nc_chunked = weighted.ComputeNC(n);

        const gauss_chunked = weighted.ComputeG(n);

        const hybrid_chunked = utils.ChunkedSum(
            hybrid_generator,
            .{ .n = n },
            0,
            n,
            utils.kChunkSize,
        );

        try writer.print("{d},,{e},\n", .{
            3 * n,
            @abs(config.kTaskAnswerP - gauss_chunked),
        });

        try writer.print("{d},{e},,\n", .{
            2 * n + 1,
            @abs(config.kTaskAnswerP - nc_chunked),
        });

        try writer.print("{d},,,{e}\n", .{
            5 * n + 1,
            @abs(config.kTaskAnswerP - hybrid_chunked),
        });

        try config.stdout.print(
            "Computed N={d:10}, Ni={d:3}/{d}\n",
            .{
                n,
                ni + 1,
                config.kSections.len,
            },
        );

        if (std.math.isNan(nc_chunked)) {
            try config.stdout.print("NewtonCotes failed\n", .{});
        }

        if (std.math.isNan(gauss_chunked)) {
            try config.stdout.print("Gauss failed\n", .{});
        }

        if (std.math.isNan(hybrid_chunked)) {
            try config.stdout.print("Hybrid failed\n", .{});
        }
    }

    try config.stdout.print("\x1B[34mWrote to output/part12.csv\x1B[0m\n", .{});
    try config.stdout.print("\x1B[1;32m<<< End\x1B[0m\n\n", .{});
}
