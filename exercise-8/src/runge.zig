const std = @import("std");
const config = @import("config.zig");

const Scalar = config.Scalar;

pub fn Norm2(x: [2]Scalar) Scalar {
    return @sqrt(x[0] * x[0] + x[1] * x[1]);
}

pub const OneStageRunge = struct {
    pub const name = "OneStageRunge";
    pub const order = 1;
    pub const calls = 1;

    pub fn call(
        _: @This(),
        f: anytype,
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
};

pub const TwoStageRunge = struct {
    pub const name = "TwoStageRunge";
    pub const order = 2;
    pub const calls = 2;

    gamma: Scalar,

    pub fn call(
        self: @This(),
        f: anytype,
        y0: [2]Scalar,
        x0: Scalar,
        step: Scalar,
    ) [2]Scalar {
        const k1 = f.call(x0, y0);
        const k2 = f.call(
            x0 + self.gamma * step,
            [2]Scalar{
                y0[0] + self.gamma * step * k1[0],
                y0[1] + self.gamma * step * k1[1],
            },
        );

        return [2]Scalar{
            y0[0] + step * ((1.0 - 1.0 / (2.0 * self.gamma)) * k1[0] +
                1.0 / (2.0 * self.gamma) * k2[0]),
            y0[1] + step * ((1.0 - 1.0 / (2.0 * self.gamma)) * k1[1] +
                1.0 / (2.0 * self.gamma) * k2[1]),
        };
    }
};

pub const ThreeStageRunge = struct {
    pub const name = "ThreeStageRunge";
    pub const order = 3;
    pub const calls = 3;

    pub fn call(
        _: @This(),
        f: anytype,
        y0: [2]Scalar,
        x0: Scalar,
        step: Scalar,
    ) [2]Scalar {
        const k1 = f.call(x0, y0);
        const k2 = f.call(
            x0 + step / 3.0,
            [2]Scalar{ y0[0] + step * k1[0] / 3.0, y0[1] + step * k1[1] / 3.0 },
        );
        const k3 = f.call(
            x0 + step * 2.0 / 3.0,
            [2]Scalar{
                y0[0] + step * k2[0] * 2.0 / 3.0,
                y0[1] + step * k2[1] * 2.0 / 3.0,
            },
        );

        return [2]Scalar{
            y0[0] + step * (k1[0] + 3.0 * k3[0]) / 4.0,
            y0[1] + step * (k1[1] + 3.0 * k3[1]) / 4.0,
        };
    }
};

pub const FourStageRunge = struct {
    pub const name = "FourStageRunge";
    pub const order = 4;
    pub const calls = 4;

    pub fn call(
        _: @This(),
        f: anytype,
        y0: [2]Scalar,
        x0: Scalar,
        step: Scalar,
    ) [2]Scalar {
        const k1 = f.call(x0, y0);
        const k2 = f.call(
            x0 + step / 2.0,
            [2]Scalar{ y0[0] + step * k1[0] / 2.0, y0[1] + step * k1[1] / 2.0 },
        );
        const k3 = f.call(
            x0 + step / 2.0,
            [2]Scalar{ y0[0] + step * k2[0] / 2.0, y0[1] + step * k2[1] / 2.0 },
        );
        const k4 = f.call(
            x0 + step,
            [2]Scalar{ y0[0] + step * k3[0], y0[1] + step * k3[1] },
        );

        return [2]Scalar{
            y0[0] + step * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0,
            y0[1] + step * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0,
        };
    }
};

fn MinStepSizeFrom(
    f: anytype,
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
    f: anytype,
    y0: [2]Scalar,
    x0: Scalar,
    order: usize,
    precision: Scalar,
) Scalar {
    const h0 = MinStepSizeFrom(f, y0, x0, order, precision);
    const y1 = (OneStageRunge{}).call(f, y0, x0, h0);
    const h1 = MinStepSizeFrom(f, y1, x0 + h0, order, precision);

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
