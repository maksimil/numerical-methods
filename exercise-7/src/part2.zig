const config = @import("config.zig");
const richardson = @import("richardson.zig");
const weighted = @import("weighted.zig");

const Scalar = config.Scalar;

pub const kTargetEps = 1e-6;

pub const NewtonCotesFormula = struct {
    pub const kConvergence = 3;

    pub fn call(n: usize) Scalar {
        return weighted.ComputeNC(n);
    }
};

pub const GaussFormula = struct {
    pub const kConvergence = 5;

    pub fn call(n: usize) Scalar {
        return weighted.ComputeG(n);
    }
};

pub fn Run() !void {
    try config.stdout.print("\x1B[1;32m>>> Part 2.1 & 2.2\x1B[0m\n", .{});

    try config.stdout.print("\x1B[34mNewton-Cotes\x1B[0m\n", .{});
    try richardson.ProcessFormula(NewtonCotesFormula, kTargetEps);

    try config.stdout.print("\x1B[34mGauss\x1B[0m\n", .{});
    try richardson.ProcessFormula(GaussFormula, kTargetEps);

    try config.stdout.print("\x1B[1;32m<<< End\x1B[0m\n\n", .{});
}
