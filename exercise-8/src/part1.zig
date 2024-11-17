const std = @import("std");
const config = @import("config.zig");
const runge = @import("runge.zig");

const Scalar = config.Scalar;

const kTaskEps = 1e-4;

fn RunMethod(
    comptime Method: type,
    method: Method,
    y0: [2]Scalar,
    x0: Scalar,
    x_end: Scalar,
    step: Scalar,
) [2]Scalar {
    var yk = y0;
    var k = @as(usize, 0);

    while (x0 + config.ToScalar(k + 1) * step < x_end) {
        yk = method.call(yk, x0 + config.ToScalar(k) * step, step);
        k += 1;
    }

    yk = method.call(
        yk,
        x0 + config.ToScalar(k) * step,
        x_end - (x0 + config.ToScalar(k) * step),
    );

    return yk;
}

pub fn Run() !void {
    try config.stdout.print("\x1B[1;32m>> Part 1\x1B[0m\n", .{});

    var step = runge.MinStepSize(
        config.TaskF,
        config.TaskF{},
        config.kTaskInitial,
        0.0,
        runge.kTwoStageRungeOrder,
        kTaskEps,
    );

    try config.stdout.print("Initial step = {e:10.3}\n", .{step});

    const run_with_step = struct {
        pub fn call(_: @This(), h: Scalar) [2]Scalar {
            const method_type = runge.TwoStageRungeMethod(config.TaskF);

            return RunMethod(
                method_type,
                method_type{ .f = config.TaskF{}, .gamma = config.kTaskXi },
                config.kTaskInitial,
                0.0,
                config.kTaskT,
                h,
            );
        }
    }{};

    var prev_result: [2]Scalar = undefined;
    var current_result: [2]Scalar = undefined;
    var loss = std.math.inf(Scalar);

    current_result = run_with_step.call(step);
    step /= 2.0;

    while (loss >= kTaskEps) {
        prev_result = current_result;
        current_result = run_with_step.call(step);

        loss = runge.RungeErr2(
            prev_result,
            current_result,
            runge.kTwoStageRungeOrder,
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
