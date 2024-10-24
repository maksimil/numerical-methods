const std = @import("std");
const config = @import("config.zig");
const weighted = @import("weighted.zig");
const part2 = @import("part2.zig");

const Scalar = config.Scalar;

const kTargetEps = 1e-6;
const kMulOrder = 2;
const kMulOrderScalar = @as(Scalar, @floatFromInt(kMulOrder));
const kLogMulOrder = @log(kMulOrderScalar);

fn Process(comptime Formula: type) !void {
    const h = config.kTaskB - config.kTaskA;

    const s1 = Formula.call(1);
    const s2 = Formula.call(kMulOrder);
    const s3 = Formula.call(kMulOrder * kMulOrder);
    const eit = -@log(@abs((s3 - s2) / (s2 - s1))) / kLogMulOrder;

    const r2 = @abs(s3 - s2) / (1 - std.math.pow(Scalar, kMulOrderScalar, -eit));

    const hopt = (h / kMulOrderScalar) * std.math.pow(Scalar, kTargetEps / r2, 1.0 / eit);
    const nopt = @as(usize, @intFromFloat(std.math.ceil(h / hopt)));

    const value = Formula.call(nopt);

    try config.stdout.print(
        "hopt={e:10.3} nopt={d:5} err={e:10.3}\n",
        .{ hopt, nopt, config.kTaskAnswerP - value },
    );
}

pub fn Run() !void {
    try config.stdout.print("\x1B[1;32m>>> Part 2.3\x1B[0m\n", .{});

    try config.stdout.print("\x1B[34mNewton-Cotes\x1B[0m\n", .{});
    try Process(part2.NewtonCotesFormula);

    try config.stdout.print("\x1B[34mGauss\x1B[0m\n", .{});
    try Process(part2.GaussFormula);

    try config.stdout.print("\x1B[1;32m<<< End\x1B[0m\n\n", .{});
}
