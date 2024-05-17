const std = @import("std");

pub fn CsvWriter() type {
    return struct {
        const Self = @This();
        const Writer = std.io.AnyWriter;

        writer_: Writer,

        pub fn init(writer: Writer) Self {
            return Self{ .writer_ = writer };
        }
    };
}
