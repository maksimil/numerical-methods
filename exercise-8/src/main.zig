const std = @import("std");
const config = @import("config.zig");
const part1 = @import("part1.zig");

pub fn main() !void {
    config.RuntimeInitialize();
    defer config.RuntimeDeinitialize();

    try part1.Run();
}
