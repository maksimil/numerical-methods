const std = @import("std");
const config = @import("config.zig");
const runge = @import("runge.zig");
const part2 = @import("part2.zig");

const Scalar = config.Scalar;

const kTaskEps = 1e-5;

const kTaskTryEps = [_]Scalar{ 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 };

fn Analyze12(method: anytype, callback: anytype) !void {
    try part2.RunMethodDynamicSteps(method, kTaskEps, struct {
        ptr_callback: @TypeOf(callback),

        pub fn call(
            self: @This(),
            y1: [2]Scalar,
            x0: Scalar,
            x1: Scalar,
            runge_err2: Scalar,
        ) !void {
            const correct = config.TaskSolution(x1);
            const err2 = runge.Norm2(
                [2]Scalar{ correct[0] - y1[0], correct[1] - y1[1] },
            );
            const err_ratio = err2 / runge_err2;

            try self.ptr_callback.call(x0, x1, err_ratio);
        }
    }{ .ptr_callback = callback });
}

pub fn Run() !void {
    try config.stdout.print("\x1B[1;32m>> Part 3.3\x1B[0m\n", .{});

    var output_dir = try std.fs.cwd().makeOpenPath("output", .{});

    const file1 = try output_dir.createFile("part331.csv", .{});
    defer file1.close();
    const writer1 = file1.writer();

    const file2 = try output_dir.createFile("part332.csv", .{});
    defer file2.close();
    const writer2 = file2.writer();

    // const file3 = try output_dir.createFile("part333.zig", .{});
    // defer file3.close();
    // const writer3 = file3.writer();

    try writer1.print("x,OneStep,TwoStep,ThreeStep,FourStep\n", .{});
    try writer2.print("x,OneRatio,TwoRatio,ThreeRatio,FourRatio\n", .{});

    try Analyze12(runge.OneStageRunge{}, struct {
        ref_writer1: @TypeOf(writer1),
        ref_writer2: @TypeOf(writer2),

        pub fn call(
            self: @This(),
            x0: Scalar,
            x1: Scalar,
            ratio: Scalar,
        ) !void {
            try self.ref_writer1.print("{e},{e},,,\n", .{ x0, x1 - x0 });
            try self.ref_writer2.print("{e},{e},,,\n", .{ x1, ratio });
        }
    }{ .ref_writer1 = writer1, .ref_writer2 = writer2 });

    try Analyze12(runge.TwoStageRunge{ .gamma = config.kTaskXi }, struct {
        ref_writer1: @TypeOf(writer1),
        ref_writer2: @TypeOf(writer2),

        pub fn call(
            self: @This(),
            x0: Scalar,
            x1: Scalar,
            ratio: Scalar,
        ) !void {
            try self.ref_writer1.print("{e},,{e},,\n", .{ x0, x1 - x0 });
            try self.ref_writer2.print("{e},,{e},,\n", .{ x1, ratio });
        }
    }{ .ref_writer1 = writer1, .ref_writer2 = writer2 });

    try Analyze12(runge.ThreeStageRunge{}, struct {
        ref_writer1: @TypeOf(writer1),
        ref_writer2: @TypeOf(writer2),

        pub fn call(
            self: @This(),
            x0: Scalar,
            x1: Scalar,
            ratio: Scalar,
        ) !void {
            try self.ref_writer1.print("{e},,,{e},\n", .{ x0, x1 - x0 });
            try self.ref_writer2.print("{e},,,{e},\n", .{ x1, ratio });
        }
    }{ .ref_writer1 = writer1, .ref_writer2 = writer2 });

    try Analyze12(runge.FourStageRunge{}, struct {
        ref_writer1: @TypeOf(writer1),
        ref_writer2: @TypeOf(writer2),

        pub fn call(
            self: @This(),
            x0: Scalar,
            x1: Scalar,
            ratio: Scalar,
        ) !void {
            try self.ref_writer1.print("{e},,,,{e}\n", .{ x0, x1 - x0 });
            try self.ref_writer2.print("{e},,,,{e}\n", .{ x1, ratio });
        }
    }{ .ref_writer1 = writer1, .ref_writer2 = writer2 });
}
