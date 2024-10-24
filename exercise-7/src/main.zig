const config = @import("config.zig");
const part11 = @import("part11.zig");
const part12 = @import("part12.zig");
const part21 = @import("part21.zig");

pub fn main() !void {
    config.RuntimeInitialize();
    defer config.RuntimeDeinitialize();

    // try part11.Run();
    // try part12.Run();
    try part21.Run();
}
