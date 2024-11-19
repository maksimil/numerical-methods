const std = @import("std");
const config = @import("config.zig");
const part1 = @import("part1.zig");
const part2 = @import("part2.zig");
const part32 = @import("part32.zig");
const part33 = @import("part33.zig");

pub fn main() !void {
    config.RuntimeInitialize();
    defer config.RuntimeDeinitialize();

    try part1.Run();
    try part2.Run();
    try part32.Run();
    try part33.Run();
}
