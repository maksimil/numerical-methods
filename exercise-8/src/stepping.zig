const std = @import("std");
const config = @import("config.zig");

const Scalar = config.Scalar;

pub fn RunSteps(
    f: anytype,
    method: anytype,
    y0: [2]Scalar,
    x0: Scalar,
    x_end: Scalar,
    step: Scalar,
) [2]Scalar {
    var yk = y0;
    var k = @as(usize, 0);

    while (x0 + config.ToScalar(k + 1) * step < x_end) {
        yk = method.call(f, yk, x0 + config.ToScalar(k) * step, step);
        k += 1;
    }

    yk = method.call(
        f,
        yk,
        x0 + config.ToScalar(k) * step,
        x_end - (x0 + config.ToScalar(k) * step),
    );

    return yk;
}
