const std = @import("std");
const config = @import("config.zig");
const matrix = @import("matrix.zig");
const lu = @import("lu.zig");
const utils = @import("utils.zig");
const least_squares = @import("least_squares.zig");

const Index = config.Index;
const Scalar = config.Scalar;
const RndGen = std.Random.DefaultPrng;

fn task_function(x: Scalar) Scalar {
    if (@abs(x) < 1e-5) {
        return @log(1e-10);
    } else {
        return @log(x * x) + (x * x * x);
    }
}

pub fn main() !void {
    // allocator
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.print("allocator: {}\n", .{gpa.deinit()});
    const allocator = gpa.allocator();

    // rng
    var rnd = RndGen.init(@bitCast(std.time.milliTimestamp()));

    // params
    const a_point: Scalar = -1;
    const b_point: Scalar = 1;
    const npoints: Index = 50;
    const data_per_point: Index = 4;
    const degrees = [_]Index{ 0, 1, 3, 5, 10 };
    const noise_magnitude: Scalar = 0.1;
    const plot_points: Index = 200;

    // creating data
    const points = try allocator.alloc(Scalar, npoints * data_per_point);
    defer allocator.free(points);

    const data = try allocator.alloc(Scalar, npoints * data_per_point);
    defer allocator.free(data);

    for (0..npoints) |i| {
        const ks = data_per_point * i;

        points[ks] = a_point +
            (b_point - a_point) *
            (@as(Scalar, @floatFromInt(i)) /
            (@as(Scalar, @floatFromInt(npoints)) - 1));

        const target_value = task_function(points[ks]);

        for (0..data_per_point) |j| {
            const k = ks + j;
            points[k] = points[ks];
            data[k] = target_value +
                utils.rand_scalar(
                rnd.random(),
                -noise_magnitude,
                noise_magnitude,
            );
        }
    }

    // computing least squares
    var coefs: [degrees.len][]Scalar = undefined;

    for (0..degrees.len) |k| {
        coefs[k] = try allocator.alloc(Scalar, degrees[k] + 1);

        const least_squares_computed = try least_squares.least_squares(
            type,
            allocator,
            coefs[k],
            utils.PolynomialCollection,
            points,
            data,
        );

        std.debug.print(
            "least_squares_computed n={}: {}\n",
            .{ degrees[k], least_squares_computed },
        );
    }

    defer {
        for (0..degrees.len) |k| {
            allocator.free(coefs[k]);
        }
    }

    // plotting
    const cwd = std.fs.cwd();
    var output_dir = try cwd.makeOpenPath("output", .{});
    defer output_dir.close();
    const graph_file = try output_dir.createFile("graph.csv", .{});
    defer graph_file.close();

    _ = try graph_file.write("x;f;data");

    for (0..degrees.len) |i| {
        const buf = try std.fmt.allocPrint(
            allocator,
            ";approx_{}",
            .{degrees[i]},
        );
        defer allocator.free(buf);
        _ = try graph_file.write(buf);
    }
    _ = try graph_file.write("\n");

    // plot the data
    for (0..npoints * data_per_point) |i| {
        const x = points[i];
        const d = data[i];

        const buf = try std.fmt.allocPrint(
            allocator,
            "{};;{}",
            .{ x, d },
        );
        defer allocator.free(buf);
        _ = try graph_file.write(buf);

        for (0..degrees.len) |_| {
            _ = try graph_file.write(";");
        }
        _ = try graph_file.write("\n");
    }

    // plot the task function
    for (0..plot_points) |i| {
        const x = a_point + (b_point - a_point) *
            (@as(Scalar, @floatFromInt(i)) /
            (@as(Scalar, @floatFromInt(plot_points)) - 1));

        const f = task_function(x);

        const buf = try std.fmt.allocPrint(
            allocator,
            "{};{};",
            .{ x, f },
        );
        defer allocator.free(buf);
        _ = try graph_file.write(buf);

        for (0..degrees.len) |k| {
            const v = utils.collection_call(
                type,
                utils.PolynomialCollection,
                coefs[k],
                x,
            );

            const buf1 = try std.fmt.allocPrint(allocator, ";{}", .{v});
            defer allocator.free(buf1);
            _ = try graph_file.write(buf1);
        }
        _ = try graph_file.write("\n");
    }
}
