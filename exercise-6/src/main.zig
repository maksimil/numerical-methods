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
    const npoints: Index = 500;
    const data_per_point: Index = 4;
    var degrees: [100]Index = undefined;
    for (0..100) |k| {
        degrees[k] = k;
    }
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

    // computing the orthogonal collection
    const OrthogonalCollection = utils.OrtogonalCollection(
        utils.InstanceCollection(utils.PolynomialCollection),
    );
    var orthogonal_collection = try OrthogonalCollection.init(
        allocator,
        degrees[degrees.len - 1] + 1,
    );
    defer orthogonal_collection.deinit(allocator);
    orthogonal_collection.compute(
        utils.InstanceCollection(utils.PolynomialCollection){},
        points,
    );

    // computing least squares
    var plain_coefs: [degrees.len][]Scalar = undefined;
    var orthogonal_coefs: [degrees.len][]Scalar = undefined;

    for (0..degrees.len) |k| {
        plain_coefs[k] = try allocator.alloc(Scalar, degrees[k] + 1);

        const CollectionType = type;
        const Orthogonal = false;
        const collection = utils.PolynomialCollection;

        // const CollectionType = utils.OrtogonalCollection(
        //     utils.InstanceCollection(utils.PolynomialCollection),
        // );
        // const Orthogonal = true;
        // var collection = try CollectionType.init(allocator, degrees[k] + 1);
        // defer collection.deinit(allocator);
        // collection.compute(
        //     utils.InstanceCollection(utils.PolynomialCollection){},
        //     points,
        // );

        const least_squares_computed = try least_squares.least_squares(
            CollectionType,
            Orthogonal,
            allocator,
            plain_coefs[k],
            collection,
            points,
            data,
        );

        std.debug.print(
            "plain least_squares_computed n={}: {}\n",
            .{ degrees[k], least_squares_computed },
        );
    }

    for (0..degrees.len) |k| {
        orthogonal_coefs[k] = try allocator.alloc(Scalar, degrees[k] + 1);

        const CollectionType = OrthogonalCollection;
        const Orthogonal = true;
        const collection = orthogonal_collection;

        const least_squares_computed = try least_squares.least_squares(
            CollectionType,
            Orthogonal,
            allocator,
            orthogonal_coefs[k],
            collection,
            points,
            data,
        );

        std.debug.print(
            "orthogonal least_squares_computed n={}: {}\n",
            .{ degrees[k], least_squares_computed },
        );
    }

    defer {
        for (0..degrees.len) |k| {
            allocator.free(plain_coefs[k]);
            allocator.free(orthogonal_coefs[k]);
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
            ";plain_{};orth_{}",
            .{ degrees[i], degrees[i] },
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
            _ = try graph_file.write(";;");
        }
        _ = try graph_file.write("\n");
    }

    // plot the task function and approximations
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
            const v_plain = utils.collection_call(
                type,
                utils.PolynomialCollection,
                plain_coefs[k],
                x,
            );

            const v_orth = utils.collection_call(
                OrthogonalCollection,
                orthogonal_collection,
                orthogonal_coefs[k],
                x,
            );

            const buf1 = try std.fmt.allocPrint(
                allocator,
                ";{};{}",
                .{ v_plain, v_orth },
            );
            defer allocator.free(buf1);
            _ = try graph_file.write(buf1);
        }
        _ = try graph_file.write("\n");
    }

    // logging the error
    const stats_file = try output_dir.createFile("stats.csv", .{});
    defer stats_file.close();

    _ = try stats_file.write("n;plain;orth\n");

    for (0..degrees.len) |i| {
        const degree = degrees[i];

        var plain_err: Scalar = 0;
        var orth_err: Scalar = 0;

        for (0..npoints * data_per_point) |j| {
            const x = points[j];
            const data_value = data[j];

            const plain_value = utils.collection_call(
                type,
                utils.PolynomialCollection,
                plain_coefs[i],
                x,
            );

            const orth_value = utils.collection_call(
                OrthogonalCollection,
                orthogonal_collection,
                orthogonal_coefs[i],
                x,
            );

            plain_err += (data_value - plain_value) * (data_value - plain_value);
            orth_err += (data_value - orth_value) * (data_value - orth_value);
        }

        const buf = try std.fmt.allocPrint(
            allocator,
            "{};{};{}\n",
            .{ degree, plain_err, orth_err },
        );
        defer allocator.free(buf);
        _ = try stats_file.write(buf);
    }
}
