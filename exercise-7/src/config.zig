const std = @import("std");

pub const Scalar = f64;
pub const kNan = std.math.nan(Scalar);

pub const kTaskA: Scalar = 1.3;
pub const kTaskB: Scalar = 2.2;
pub const kTaskAlpha: Scalar = 0.0;
pub const kTaskBeta: Scalar = 5.0 / 6.0;
pub const kTaskAnswer: Scalar = 3.0765665777241504;
pub const kTaskAnswerP: Scalar = 10.8395451094690891;

pub fn TaskF(x: Scalar) Scalar {
    return (4.0 * @cos(0.5 * x) * @exp(-5.0 * x / 4.0)) +
        (2.0 * @sin(4.5 * x) * @exp(x / 8.0)) + 2.0;
}

pub const kSections = [_]usize{
    1,
    2,
    4,

    8,
    16,
    32,
    64,
    128,
    256,
    512,
    1024,
    2048,
    4096,
    8192,
    16384,
    32768,
    65536,
    131072,
    262144,
    524288,

    // 1048576,
    // 2097152,
    // 4194304,
    // 8388608,
    // 16777216,
    // 33554432,
    // 67108864,
    // 134217728,
    // 268435456,
    // 536870912,
};

pub const MIN_PIVOT: Scalar = 1e-14;

// runtime
pub var stdout: std.fs.File.Writer = undefined;
pub var stderr: std.fs.File.Writer = undefined;
var gpa: std.heap.GeneralPurposeAllocator(.{}) = undefined;
pub var allocator: std.mem.Allocator = undefined;

pub fn RuntimeInitialize() void {
    stdout = std.io.getStdOut().writer();
    stderr = std.io.getStdErr().writer();

    gpa = @TypeOf(gpa){};
    allocator = gpa.allocator();
}

pub fn RuntimeDeinitialize() void {
    stderr.print("allocator: {}\n", .{gpa.deinit()}) catch unreachable;
}
