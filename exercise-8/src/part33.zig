const std = @import("std");
const config = @import("config.zig");
const runge = @import("runge.zig");
const part2 = @import("part2.zig");

const Scalar = config.Scalar;

const kTaskEps = 1e-5;

const kTaskTryEps = [_]Scalar{ 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7 };

fn Analyze12(method: anytype, callback: anytype) !void {
    var callback_wrapper = struct {
        ptr_callback: @TypeOf(callback),
        y0: [2]Scalar,

        pub fn call(
            self: *@This(),
            y1: [2]Scalar,
            x0: Scalar,
            x1: Scalar,
            runge_err2: Scalar,
        ) !void {
            const correct = config.TaskSolutionIC(self.y0, x0, x1);
            const err2 = runge.Norm2(
                [2]Scalar{ correct[0] - y1[0], correct[1] - y1[1] },
            );
            const err_ratio = err2 / runge_err2;

            try self.ptr_callback.call(x0, x1, err_ratio);

            self.y0 = y1;
        }
    }{ .ptr_callback = callback, .y0 = config.kTaskInitial };

    try part2.RunMethodDynamicSteps(method, kTaskEps, &callback_wrapper);
}

fn Run12() !void {
    try config.stdout.print("\x1B[1;32m>> Part 3.3.1-2\x1B[0m\n", .{});

    var output_dir = try std.fs.cwd().makeOpenPath("output", .{});

    const file1 = try output_dir.createFile("part331.csv", .{});
    defer file1.close();
    const writer1 = file1.writer();

    const file2 = try output_dir.createFile("part332.csv", .{});
    defer file2.close();
    const writer2 = file2.writer();

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

fn CountCalls(method: anytype, precision: Scalar) !usize {
    var ncalls = @as(usize, 0);

    const method_wrapper = struct {
        pub const name = @TypeOf(method).name;
        pub const order = @TypeOf(method).order;

        ptr_method: *const @TypeOf(method),
        ptr_ncalls: *usize,

        pub fn call(
            self: @This(),
            f: anytype,
            y0: [2]Scalar,
            x0: Scalar,
            step: Scalar,
        ) [2]Scalar {
            self.ptr_ncalls.* += @as(usize, @TypeOf(method).calls);
            return self.ptr_method.call(f, y0, x0, step);
        }
    }{ .ptr_ncalls = &ncalls, .ptr_method = &method };

    try part2.RunMethodDynamicSteps(method_wrapper, precision, part2.kLoggerCallback);

    return ncalls;
}

fn Run3() !void {
    try config.stdout.print("\x1B[1;32m>> Part 3.3.1-2\x1B[0m\n", .{});

    var output_dir = try std.fs.cwd().makeOpenPath("output", .{});

    const file3 = try output_dir.createFile("part333.csv", .{});
    defer file3.close();
    const writer3 = file3.writer();

    try writer3.print("Precision,OneCalls,TwoCalls,ThreeCalls,FourCalls\n", .{});

    for (0..kTaskTryEps.len) |k| {
        const precision = kTaskTryEps[k];

        const one_calls = try CountCalls(runge.OneStageRunge{}, precision);
        const two_calls = try CountCalls(
            runge.TwoStageRunge{ .gamma = config.kTaskXi },
            precision,
        );
        const three_calls = try CountCalls(runge.ThreeStageRunge{}, precision);
        const four_calls = try CountCalls(runge.FourStageRunge{}, precision);

        try writer3.print(
            "{e},{d},{d},{d},{d}\n",
            .{ precision, one_calls, two_calls, three_calls, four_calls },
        );
    }
}

pub fn Run() !void {
    try Run12();
    try Run3();
}
