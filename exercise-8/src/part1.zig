const std = @import("std");
const config = @import("config.zig");
const runge = @import("runge.zig");
const stepping = @import("stepping.zig");

const Scalar = config.Scalar;

const kTaskEps = 1e-4;

fn RunMethod(method: anytype) !void {
    try config.stdout.print("\x1B[34mUsing {s}\x1B[0m\n", .{@TypeOf(method).name});

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

    while (loss >= kTaskEps) {
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

        step /= 2.0;
    }

    const correct_value = config.TaskSolution(config.kTaskT);

    try config.stdout.print(
        "{s:13}ans_y=[{e:10.3}, {e:10.3}]\n",
        .{ "", correct_value[0], correct_value[1] },
    );

    try config.stdout.print("\n", .{});
}

pub fn Run() !void {
    try config.stdout.print("\x1B[1;32m>> Part 1\x1B[0m\n", .{});

    try RunMethod(runge.OneStageRunge{});
    try RunMethod(runge.TwoStageRunge{ .gamma = config.kTaskXi });
    try RunMethod(runge.ThreeStageRunge{});
    try RunMethod(runge.FourStageRunge{});
}
