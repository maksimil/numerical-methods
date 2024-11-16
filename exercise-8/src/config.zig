const std = @import("std");

pub const Scalar = f64;
pub const kNan = std.math.nan(Scalar);

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

pub fn ToScalar(x: anytype) Scalar {
    return @as(Scalar, @floatFromInt(x));
}

// task
pub const kTaskXi = 1.0 / 15.0;
pub const kTaskA = 1.0 / 10.0;
pub const kTaskB = 1.0 / 12.0;

pub const kTaskT = std.math.pi;
pub const kTaskInitial = [2]Scalar{ kTaskB * std.math.pi, kTaskA * std.math.pi };

pub const kOmega = @sqrt(kTaskA * kTaskB);
pub const kLambda = @sqrt(kTaskA / kTaskB);

pub const TaskF = struct {
    pub fn call(_: @This(), _: Scalar, y: [2]Scalar) [2]Scalar {
        return [2]Scalar{ kTaskA * y[1], -kTaskB * y[0] };
    }
};

pub fn TaskSolution(x: Scalar) [2]Scalar {
    return [2]Scalar{
        std.math.pi * (kTaskB * @cos(kOmega * x) +
            kTaskA * kLambda * @sin(kOmega * x)),
        std.math.pi * (-kTaskB / kLambda * @sin(kOmega * x) +
            kTaskA * @cos(kOmega * x)),
    };
}
