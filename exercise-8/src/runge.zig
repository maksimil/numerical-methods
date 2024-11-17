const std = @import("std");
const config = @import("config.zig");

const Scalar = config.Scalar;

pub fn Norm2(x: [2]Scalar) Scalar {
    return @sqrt(x[0] * x[0] + x[1] * x[1]);
}

pub fn OneStageRunge(
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

pub fn OneStageRungeMethod(Function: type) type {
    return struct {
        f: Function,

        pub fn call(self: @This(), y0: [2]Scalar, x0: Scalar, h: Scalar) [2]Scalar {
            return OneStageRunge(Function, self.f, y0, x0, h);
        }
    };
}

pub fn TwoStageRunge(
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

pub fn TwoStageRungeMethod(Function: type) type {
    return struct {
        f: Function,
        gamma: Scalar,

        pub fn call(self: @This(), y0: [2]Scalar, x0: Scalar, h: Scalar) [2]Scalar {
            return TwoStageRunge(Function, self.f, y0, x0, h, self.gamma);
        }
    };
}

pub const kTwoStageRungeOrder = 2;

fn MinStepSizeFrom(
    comptime Function: type,
    f: Function,
    y0: [2]Scalar,
    x0: Scalar,
    order: usize,
    precision: Scalar,
) Scalar {
    const initial_f = f.call(x0, y0);

    const delta = 1.0 / std.math.pow(
        Scalar,
        config.kTaskT,
        config.ToScalar(order + 1),
    ) +
        std.math.pow(
        Scalar,
        Norm2(initial_f),
        config.ToScalar(order + 1),
    );

    const initial_step =
        std.math.pow(Scalar, precision / delta, 1.0 / config.ToScalar(order + 1));

    return initial_step;
}

pub fn MinStepSize(
    comptime Function: type,
    f: Function,
    y0: [2]Scalar,
    x0: Scalar,
    order: usize,
    precision: Scalar,
) Scalar {
    const h0 = MinStepSizeFrom(Function, f, y0, x0, order, precision);
    const y1 = OneStageRunge(Function, f, y0, x0, h0);
    const h1 = MinStepSizeFrom(Function, f, y1, x0 + h0, order, precision);

    return @min(h0, h1);
}

pub fn RungeErr2(
    prev_result: [2]Scalar,
    current_result: [2]Scalar,
    order: usize,
    ratio: Scalar,
) Scalar {
    const denominator = std.math.pow(Scalar, ratio, config.ToScalar(order)) - 1.0;
    const runge_err = [2]Scalar{
        (current_result[0] - prev_result[0]) / denominator,
        (current_result[1] - prev_result[1]) / denominator,
    };
    return Norm2(runge_err);
}
