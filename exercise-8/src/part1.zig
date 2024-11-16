const std = @import("std");
const config = @import("config.zig");

const Scalar = config.Scalar;

const kTaskEps = 1e-4;

const kTwoStageRungeOrder = 2.0;

fn OneStageRunge(
    comptime Function: type,
    f: Function,
    y0: [2]Scalar,
    x0: Scalar,
    step: Scalar,
) [2]Scalar {
    const k1 = f.call(x0, y0);

    var y1: [2]Scalar = undefined;
    y1[0] = y0[0] + step * k1[0];
    y1[1] = y0[1] + step * k1[1];

    return y1;
}

fn TwoStageRunge(
    comptime Function: type,
    f: Function,
    y0: [2]Scalar,
    x0: Scalar,
    step: Scalar,
    gamma: Scalar,
) [2]Scalar {
    const k1 = f.call(x0, y0);
    const k2 = f.call(
        x0 + gamma * step,
        [2]Scalar{ y0[0] + gamma * step * k1[0], y0[1] + gamma * step * k1[1] },
    );

    return [2]Scalar{
        y0[0] + step * ((1.0 - 1.0 / (2.0 * gamma)) * k1[0] +
            1.0 / (2.0 * gamma) * k2[0]),
        y0[1] + step * ((1.0 - 1.0 / (2.0 * gamma)) * k1[1] +
            1.0 / (2.0 * gamma) * k2[1]),
    };
}

fn TwoStageRungeMethod(Function: type) type {
    return struct {
        f: Function,
        gamma: Scalar,

        pub fn call(self: @This(), y0: [2]Scalar, x0: Scalar, h: Scalar) [2]Scalar {
            return TwoStageRunge(Function, self.f, y0, x0, h, self.gamma);
        }
    };
}

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

        // const v = config.TaskSolution(x0 + config.ToScalar(k) * step);
        // config.stdout.print(
        //     "t={e:10.3}, y=[{e:10.3}, {e:10.3}]\n" ++
        //         "        , y_ans=[{e:10.3}, {e:10.3}]\n",
        //     .{ x0 + config.ToScalar(k) * step, yk[0], yk[1], v[0], v[1] },
        // ) catch unreachable;
    }

    yk = method.call(
        yk,
        x0 + config.ToScalar(k) * step,
        x_end - (x0 + config.ToScalar(k) * step),
    );

    // const v = config.TaskSolution(x_end);
    // config.stdout.print(
    //     "t={e:10.3}, y=[{e:10.3}, {e:10.3}]\n" ++
    //         "        , y_ans=[{e:10.3}, {e:10.3}]\n",
    //     .{ x_end, yk[0], yk[1], v[0], v[1] },
    // ) catch unreachable;

    return yk;
}

fn MinStepSizeFrom(
    comptime Function: type,
    f: Function,
    y0: [2]Scalar,
    x0: Scalar,
) Scalar {
    const initial_f = f.call(x0, y0);

    const delta = 1.0 / std.math.pow(Scalar, config.kTaskT, kTwoStageRungeOrder + 1.0) +
        std.math.pow(
        Scalar,
        @sqrt(initial_f[0] * initial_f[0] + initial_f[1] * initial_f[1]),
        kTwoStageRungeOrder + 1.0,
    );

    const initial_step =
        std.math.pow(Scalar, kTaskEps / delta, 1.0 / (kTwoStageRungeOrder + 1.0));

    return initial_step;
}

fn MinStepSize(
    comptime Function: type,
    f: Function,
    y0: [2]Scalar,
    x0: Scalar,
) Scalar {
    const h0 = MinStepSizeFrom(Function, f, y0, x0);
    const y1 = OneStageRunge(Function, f, y0, x0, h0);
    const h1 = MinStepSizeFrom(Function, f, y1, x0 + h0);

    return @min(h0, h1);
}

pub fn Run() !void {
    try config.stdout.print("\x1B[1;32m>> Part 1\x1B[0m\n", .{});

    const initial_step = MinStepSize(
        config.TaskF,
        config.TaskF{},
        config.kTaskInitial,
        0.0,
    );
    var step = initial_step;

    try config.stdout.print("Initial step = {e:10.3}\n", .{initial_step});

    const run_with_step = struct {
        pub fn call(_: @This(), h: Scalar) [2]Scalar {
            const method_type = TwoStageRungeMethod(config.TaskF);

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

        const runge_err = [2]Scalar{
            (current_result[0] - prev_result[0]) /
                (std.math.pow(Scalar, 2.0, kTwoStageRungeOrder) - 1),
            (current_result[1] - prev_result[1]) /
                (std.math.pow(Scalar, 2.0, kTwoStageRungeOrder) - 1),
        };

        loss = @sqrt(runge_err[0] * runge_err[0] + runge_err[1] * runge_err[1]);
        const v = config.TaskSolution(config.kTaskT);
        const err2 = @sqrt((v[0] - current_result[0]) * (v[0] - current_result[0]) +
            (v[1] - current_result[1]) * (v[1] - current_result[1]));

        try config.stdout.print(
            "step={e:10.3}, y=[{e:10.3}, {e:10.3}], " ++
                "runge_err2={e:10.3}, err2={e:10.3}\n",
            .{ step, current_result[0], current_result[1], loss, err2 },
        );

        step /= 2.0;
    }

    const correct_value = config.TaskSolution(config.kTaskT);

    try config.stdout.print(
        "{s:15}, y=[{e:10.3}, {e:10.3}]\n",
        .{ "", correct_value[0], correct_value[1] },
    );
}
