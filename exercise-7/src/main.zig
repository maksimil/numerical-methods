const config = @import("config.zig");
const part11 = @import("part11.zig");
const part12 = @import("part12.zig");
const part2 = @import("part2.zig");
const part23 = @import("part23.zig");

pub fn main() !void {
    config.RuntimeInitialize();
    defer config.RuntimeDeinitialize();

    try part11.Run();
    try part12.Run();
    try part2.Run();
    try part23.Run();
}
