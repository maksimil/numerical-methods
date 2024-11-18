const std = @import("std");
const config = @import("config.zig");
const runge = @import("runge.zig");

const Scalar = config.Scalar;

const kTaskEps = 1e-5;

fn RunMethod(method: anytype) !void {
    try config.stdout.print("\x1B[34mUsing {s}\x1B[0m\n", .{@TypeOf(method).name});

    var step = runge.MinStepSize(
        config.kTaskF,
        config.kTaskInitial,
        0.0,
        @TypeOf(method).order,
        kTaskEps,
    );

    const err_ratio = std.math.pow(
        Scalar,
        2.0,
        -config.ToScalar(@TypeOf(method).order),
    );

    var current_y = config.kTaskInitial;
    var current_x = @as(Scalar, 0.0);

    var continue_iterations = true;
    while (continue_iterations) {
        if (config.kTaskT - current_x < step) {
            step = config.kTaskT - current_x;
            continue_iterations = false;
        }

        var runge_err2 = @as(Scalar, 0.0);
        var chosen_step = @as(Scalar, 0.0);
        while (true) {
            const y_single = method.call(
                config.kTaskF,
                current_y,
                current_x,
                step,
            );
            const y_half_first = method.call(
                config.kTaskF,
                current_y,
                current_x,
                step / 2.0,
            );
            const y_half_second = method.call(
                config.kTaskF,
                y_half_first,
                current_x + step / 2.0,
                step / 2.0,
            );

            runge_err2 = runge.RungeErr2(
                y_single,
                y_half_second,
                @TypeOf(method).order,
                2.0,
            );

            if (runge_err2 > kTaskEps) {
                step /= 2.0;
                continue_iterations = true;
            } else if (runge_err2 > kTaskEps * err_ratio) {
                current_y = y_half_second;
                current_x += step;
                chosen_step = step;
                step /= 2.0;
                break;
            } else if (runge_err2 > kTaskEps * err_ratio * err_ratio) {
                current_y = y_single;
                current_x += step;
                chosen_step = step;
                break;
            } else {
                current_y = y_single;
                current_x += step;
                chosen_step = step;
                step *= 2.0;
                break;
            }
        }

        const answer_y = config.TaskSolution(current_x);
        const err2 = runge.Norm2(
            [2]Scalar{ current_y[0] - answer_y[0], current_y[1] - answer_y[1] },
        );
        try config.stdout.print(
            "x={e:10.3}, y=[{e:10.3}, {e:10.3}], err2={e:10.3}, " ++
                "runge_err2_local={e:10.3}, step={e:10.3}\n",
            .{ current_x, current_y[0], current_y[1], err2, runge_err2, chosen_step },
        );
    }

    const answer_y = config.TaskSolution(config.kTaskT);
    const err2 = runge.Norm2(
        [2]Scalar{ current_y[0] - answer_y[0], current_y[1] - answer_y[1] },
    );

    try config.stdout.print(
        "y    =[{e:10.3}, {e:10.3}], err2={e:10.3}\n" ++
            "ans_y=[{e:10.3}, {e:10.3}]\n",
        .{ current_y[0], current_y[1], err2, answer_y[0], answer_y[1] },
    );

    try config.stdout.print("\n", .{});
}

pub fn Run() !void {
    try config.stdout.print("\x1B[1;32m>> Part 2\x1B[0m\n", .{});

    try RunMethod(runge.OneStageRunge{});
    try RunMethod(runge.TwoStageRunge{ .gamma = config.kTaskXi });
    try RunMethod(runge.ThreeStageRunge{});
    try RunMethod(runge.FourStageRunge{});
}
