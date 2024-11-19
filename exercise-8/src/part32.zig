const std = @import("std");
const config = @import("config.zig");
const runge = @import("runge.zig");
const stepping = @import("stepping.zig");

const Scalar = config.Scalar;

const kTaskEps = 1e-4;

const kTrueValueStep = 1e-4;

fn FindMethodStep(method: anytype) !Scalar {
    var step = runge.MinStepSize(
        config.kTaskF,
        config.kTaskInitial,
        0.0,
        @TypeOf(method).order,
        kTaskEps,
    );

    try config.stdout.print("Initial step = {e:10.3}\n", .{step});

    var prev_result: [2]Scalar = undefined;
    var current_result: [2]Scalar = undefined;
    var loss = std.math.inf(Scalar);

    current_result = stepping.RunSteps(
        config.kTaskF,
        method,
        config.kTaskInitial,
        0.0,
        config.kTaskT,
        step,
    );
    step /= 2.0;

    while (true) {
        prev_result = current_result;
        current_result = stepping.RunSteps(
            config.kTaskF,
            method,
            config.kTaskInitial,
            0.0,
            config.kTaskT,
            step,
        );

        loss = runge.RungeErr2(
            prev_result,
            current_result,
            @TypeOf(method).order,
            2.0,
        );
        const v = config.TaskSolution(config.kTaskT);
        const err2 = runge.Norm2(
            [2]Scalar{ v[0] - current_result[0], v[1] - current_result[1] },
        );

        try config.stdout.print(
            "step={e:10.3}, y=[{e:10.3}, {e:10.3}], " ++
                "runge_err2={e:10.3}, err2={e:10.3}\n",
            .{ step, current_result[0], current_result[1], loss, err2 },
        );

        if (loss < kTaskEps) {
            break;
        } else {
            step /= 2.0;
        }
    }

    return step;
}

fn LogPoints(method: anytype, step: Scalar, log_point: anytype) !void {
    const y0 = config.kTaskInitial;
    const x0 = 0.0;
    const x_end = config.kTaskT;
    const f = config.kTaskF;

    var yk = y0;
    var k = @as(usize, 0);

    try log_point.call(y0, x0);

    while (x0 + config.ToScalar(k + 1) * step < x_end) {
        yk = method.call(f, yk, x0 + config.ToScalar(k) * step, step);
        k += 1;

        try log_point.call(yk, x0 + config.ToScalar(k) * step);
    }

    yk = method.call(
        f,
        yk,
        x0 + config.ToScalar(k) * step,
        x_end - (x0 + config.ToScalar(k) * step),
    );

    try log_point.call(yk, x_end);
}

fn RunMethod(method: anytype, log_point: anytype) !void {
    try config.stdout.print("\x1B[34mUsing {s}\x1B[0m\n", .{@TypeOf(method).name});

    const step = try FindMethodStep(method);
    try LogPoints(method, step, log_point);

    try config.stdout.print("\n", .{});
}

fn Err2(y: [2]Scalar, x: Scalar) Scalar {
    const z = config.TaskSolution(x);

    return runge.Norm2([2]Scalar{ y[0] - z[0], y[1] - z[1] });
}

pub fn Run() !void {
    try config.stdout.print("\x1B[1;32m>> Part 3.2\x1B[0m\n", .{});

    var output_dir = try std.fs.cwd().makeOpenPath("output", .{});
    const file = try output_dir.createFile("part32.csv", .{});
    defer file.close();
    const writer = file.writer();

    try writer.print("x,OneErr,TwoErr,ThreeErr,FourErr\n", .{});

    try RunMethod(runge.OneStageRunge{}, struct {
        ptr_writer: *const @TypeOf(writer),

        pub fn call(self: @This(), y: [2]Scalar, x: Scalar) !void {
            try self.ptr_writer.print("{e},{e},,,\n", .{ x, Err2(y, x) });
        }
    }{ .ptr_writer = &writer });

    try RunMethod(runge.TwoStageRunge{ .gamma = config.kTaskXi }, struct {
        ptr_writer: *const @TypeOf(writer),

        pub fn call(self: @This(), y: [2]Scalar, x: Scalar) !void {
            try self.ptr_writer.print("{e},,{e},,\n", .{ x, Err2(y, x) });
        }
    }{ .ptr_writer = &writer });

    try RunMethod(runge.ThreeStageRunge{}, struct {
        ptr_writer: *const @TypeOf(writer),

        pub fn call(self: @This(), y: [2]Scalar, x: Scalar) !void {
            try self.ptr_writer.print("{e},,,{e},\n", .{ x, Err2(y, x) });
        }
    }{ .ptr_writer = &writer });

    try RunMethod(runge.FourStageRunge{}, struct {
        ptr_writer: *const @TypeOf(writer),

        pub fn call(self: @This(), y: [2]Scalar, x: Scalar) !void {
            try self.ptr_writer.print("{e},,,,{e}\n", .{ x, Err2(y, x) });
        }
    }{ .ptr_writer = &writer });
}
