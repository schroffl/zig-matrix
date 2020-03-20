const std = @import("std");
const testing = std.testing;

pub const StorageOrder = enum {
    RowMajor,
    ColumnMajor,
};

fn MultiplicationType(comptime A: type, comptime B: type) type {
    comptime std.debug.assert(A.Child == B.Child);
    comptime if (A.columns != B.rows) {
        @compileError("MxN matrices can only be multiplied with NxP matrices");
    };

    return Matrix(A.Child, A.rows, B.columns);
}

fn TranspositionType(comptime T: type) type {
    return Matrix(T.Child, T.columns, T.rows);
}

pub fn Matrix(
    comptime T: type,
    comptime m: usize,
    comptime n: usize,
) type {
    return struct {
        const Self = @This();

        pub const Child = T;
        pub const rows = m;
        pub const columns = n;

        data: [rows * columns]T,

        pub fn init() Self {
            return Self{
                .data = [_]T{0} ** (rows * columns),
            };
        }

        pub fn fromValues(data: [rows * columns]T, order: StorageOrder) Self {
            var out = Self.init();
            out.assign(data, order);
            return out;
        }

        pub usingnamespace if (comptime rows == columns)
            struct {
                pub const identity = comptime blk: {
                    var mat = comptime Self.init();
                    var i = 0;

                    while (i < rows) : (i += 1) {
                        mat.set(i, i, 1);
                    }

                    break :blk mat;
                };

                pub fn determinant(self: Self) T {
                    switch (rows) {
                        2 => {
                            return self.get(0, 0) * self.get(1, 1) - self.get(0, 1) * self.get(1, 0);
                        },
                        3 => {
                            const a00 = self.get(0, 0);
                            const a01 = self.get(0, 1);
                            const a02 = self.get(0, 2);

                            const a10 = self.get(1, 0);
                            const a11 = self.get(1, 1);
                            const a12 = self.get(1, 2);

                            const a20 = self.get(2, 0);
                            const a21 = self.get(2, 1);
                            const a22 = self.get(2, 2);

                            return (a00 * a11 * a22) + (a01 * a12 * a20) + (a02 * a10 * a21) - (a20 * a11 * a02) - (a21 * a12 * a00) - (a22 * a10 * a01);
                        },
                        else => @compileError("TODO Implement for types other than 2x2 3x3"),
                    }
                }

                pub fn invert(self: Self) error{SingularMatrix}!Self {
                    const det = self.determinant();

                    if (det == 0) {
                        return error.SingularMatrix;
                    }

                    var out = Self.init();

                    switch (rows) {
                        2 => {
                            const a00 = self.get(0, 0);
                            const a11 = self.get(1, 1);

                            const a01 = self.get(0, 1);
                            const a10 = self.get(1, 0);

                            out.set(0, 0, a11 * det);
                            out.set(1, 1, a00 * det);

                            out.set(0, 1, -a01 * det);
                            out.set(1, 0, -a10 * det);
                        },
                        3 => {
                            const a00 = self.get(0, 0);
                            const a01 = self.get(0, 1);
                            const a02 = self.get(0, 2);

                            const a10 = self.get(1, 0);
                            const a11 = self.get(1, 1);
                            const a12 = self.get(1, 2);

                            const a20 = self.get(2, 0);
                            const a21 = self.get(2, 1);
                            const a22 = self.get(2, 2);

                            out.set(0, 0, a11 * a22 - a12 * a21);
                            out.set(0, 1, a02 * a21 - a01 * a22);
                            out.set(0, 2, a01 * a12 - a02 * a11);

                            out.set(1, 0, a12 * a20 - a10 * a22);
                            out.set(1, 1, a00 * a22 - a02 * a20);
                            out.set(1, 2, a02 * a10 - a00 * a12);

                            out.set(2, 0, a10 * a21 - a11 * a20);
                            out.set(2, 1, a01 * a20 - a00 * a21);
                            out.set(2, 2, a00 * a11 - a01 * a10);

                            out = out.multiplyScalar(1 / det);
                        },
                        else => @compileError("TODO Implement for types other than 2x2 3x3"),
                    }

                    return out;
                }
            }
        else
        // Non-quadratic matrices don't have an identity.
            struct {};

        pub fn add(a: Self, b: Self) Self {
            var out = Self.init();
            var i: usize = 0;

            while (i < rows) : (i += 1) {
                var j: usize = 0;

                while (j < columns) : (j += 1) {
                    const value = a.get(i, j) + b.get(i, j);
                    out.set(i, j, value);
                }
            }

            return out;
        }

        pub fn subtract(a: Self, b: Self) Self {
            var out = Self.init();
            var i: usize = 0;

            while (i < rows) : (i += 1) {
                var j: usize = 0;

                while (j < columns) : (j += 1) {
                    out.set(i, j, a.get(i, j) - b.get(i, j));
                }
            }

            return out;
        }

        pub fn multiplyScalar(a: Self, scalar: T) Self {
            var out = Self.init();
            var i: usize = 0;

            while (i < rows) : (i += 1) {
                var j: usize = 0;

                while (j < columns) : (j += 1) {
                    out.set(i, j, a.get(i, j) * scalar);
                }
            }

            return out;
        }

        pub fn multiply(a: Self, b: var) MultiplicationType(@TypeOf(a), @TypeOf(b)) {
            comptime var Out = MultiplicationType(@TypeOf(a), @TypeOf(b));
            var out = Out.init();

            var i: usize = 0;

            while (i < Out.rows) : (i += 1) {
                var j: usize = 0;

                while (j < Out.columns) : (j += 1) {
                    var k: usize = 0;
                    var result: Out.Child = 0;

                    while (k < @TypeOf(a).columns) {
                        result += a.get(i, k) * b.get(k, j);
                        k += 1;
                    }

                    out.set(i, j, result);
                }
            }

            return out;
        }

        pub inline fn get(self: Self, row: usize, column: usize) T {
            std.debug.assert(row < rows);
            std.debug.assert(column < columns);

            return self.data[row * columns + column];
        }

        pub inline fn set(self: *Self, row: usize, column: usize, value: T) void {
            std.debug.assert(row < rows);
            std.debug.assert(column < columns);

            self.data[row * columns + column] = value;
        }

        pub fn assign(self: *Self, data: [rows * columns]T, order: StorageOrder) void {
            switch (order) {
                .RowMajor => std.mem.copy(T, &self.data, &data),
                .ColumnMajor => {
                    var i: usize = 0;

                    while (i < rows) : (i += 1) {
                        var j: usize = 0;

                        while (j < columns) : (j += 1) {
                            self.set(i, j, data[j * rows + i]);
                        }
                    }
                },
            }
        }

        pub fn getData(self: Self, order: StorageOrder) [rows * columns]T {
            switch (order) {
                .RowMajor => return self.data,
                .ColumnMajor => {
                    var data = [_]T{0} ** (rows * columns);

                    var i: usize = 0;

                    while (i < rows) : (i += 1) {
                        var j: usize = 0;

                        while (j < columns) : (j += 1) {
                            data[j * rows + i] = self.get(i, j);
                        }
                    }

                    return data;
                },
            }
        }

        pub fn transpose(self: Self) TranspositionType(Self) {
            var out = TranspositionType(Self).init();
            var i: usize = 0;

            while (i < rows) : (i += 1) {
                var j: usize = 0;

                while (j < columns) : (j += 1) {
                    out.set(j, i, self.get(i, j));
                }
            }

            return out;
        }

        pub fn eql(a: Self, b: var) bool {
            var i: usize = 0;

            comptime std.debug.assert(@TypeOf(b).columns == Self.columns);
            comptime std.debug.assert(@TypeOf(b).rows == Self.rows);

            while (i < m) : (i += 1) {
                var j: usize = 0;

                while (j < n) : (j += 1) {
                    if (a.get(i, j) != b.get(i, j)) {
                        return false;
                    }
                }
            }

            return true;
        }

        pub fn approxEql(a: Self, b: var, epsilon: T) bool {
            var i: usize = 0;

            comptime std.debug.assert(@TypeOf(b).columns == Self.columns);
            comptime std.debug.assert(@TypeOf(b).rows == Self.rows);

            while (i < m) : (i += 1) {
                var j: usize = 0;

                while (j < n) : (j += 1) {
                    if (!std.math.approxEq(T, a.get(i, j), b.get(i, j), epsilon)) {
                        return false;
                    }
                }
            }

            return true;
        }

        pub fn format(
            self: Self,
            comptime fmt: []const u8,
            options: std.fmt.FormatOptions,
            out_stream: var,
        ) @TypeOf(out_stream).Error!void {
            try out_stream.writeAll("Matrix<");
            try std.fmt.format(out_stream, "{}", .{rows});
            try out_stream.writeAll("x");
            try std.fmt.format(out_stream, "{}", .{columns});
            try out_stream.writeAll(">");
            try out_stream.writeAll("(");
            try std.fmt.format(out_stream, "{}", .{@typeName(T)});
            try out_stream.writeAll(")\n");

            var i: usize = 0;

            while (i < Self.rows) : (i += 1) {
                var j: usize = 0;

                // TODO Align the columns and print column indices

                if (i == 0) {
                    try std.fmt.format(out_stream, "{} \u{23A1} ", .{i});
                } else if (i == rows - 1) {
                    try std.fmt.format(out_stream, "{} \u{23A3} ", .{i});
                } else {
                    try std.fmt.format(out_stream, "{} \u{23A2} ", .{i});
                }

                while (j < Self.columns) : (j += 1) {
                    if (j > 0) {
                        try out_stream.writeAll(", ");
                    }

                    try std.fmt.formatType(self.get(i, j), fmt, options, out_stream, 1);
                }

                if (i == 0) {
                    try out_stream.writeAll(" \u{23A4}\n");
                } else if (i == rows - 1) {
                    try out_stream.writeAll(" \u{23A6}\n");
                } else {
                    try out_stream.writeAll(" \u{23A5}\n");
                }
            }
        }
    };
}

test "Equality" {
    const Mat2 = Matrix(u8, 2, 2);

    const a = Mat2.fromValues([_]u8{
        5,  34,
        20, 4,
    }, .RowMajor);

    const b = Mat2.fromValues([_]u8{
        5,  34,
        20, 4,
    }, .RowMajor);

    const c = Mat2.fromValues([_]u8{
        4,  34,
        20, 5,
    }, .RowMajor);

    testing.expect(a.eql(b));
    testing.expect(b.eql(a));
    testing.expect(!a.eql(c));
    testing.expect(!b.eql(c));
    testing.expect(!c.eql(a));
    testing.expect(!c.eql(b));

    const Mat2Float = Matrix(f32, 2, 2);
    const d = Mat2Float.fromValues([_]f32{
        5,  34.134234,
        20, 4,
    }, .RowMajor);

    testing.expect(d.eql(d));
}

test "Approximate Equality" {
    const Mat2 = Matrix(f32, 2, 2);

    const a = Mat2.fromValues([_]f32{
        5,  34.134234,
        20, 4,
    }, .RowMajor);

    const b = Mat2.fromValues([_]f32{
        5,  34.134231,
        20, 4,
    }, .RowMajor);

    testing.expect(a.approxEql(b, 0.000004));
    testing.expect(!a.approxEql(b, 0.000003));

    testing.expect(b.approxEql(a, 0.000004));
    testing.expect(!b.approxEql(a, 0.000003));
}

test "Addition / Subtraction" {
    const Mat2 = Matrix(f32, 2, 2);
    const id = Mat2.identity;
    const mat = Mat2.init();

    const addResult = Mat2.add(mat, id);
    testing.expect(id.eql(addResult));

    const subResult = Mat2.subtract(mat, id);
    testing.expect(id.multiplyScalar(-1).eql(subResult));
}

test "Transposing a 2x4 Matrix" {
    const Mat2x4 = Matrix(f32, 2, 4);
    var mat = Mat2x4.init();

    mat.set(0, 0, 5);
    mat.set(0, 3, 10);

    var transposed = mat.transpose();

    testing.expectEqual(@as(usize, 2), Mat2x4.rows);
    testing.expectEqual(@as(usize, 4), Mat2x4.columns);

    testing.expectEqual(@as(usize, 4), @TypeOf(transposed).rows);
    testing.expectEqual(@as(usize, 2), @TypeOf(transposed).columns);

    testing.expectEqual(mat.get(0, 0), transposed.get(0, 0));
    testing.expectEqual(mat.get(0, 3), transposed.get(3, 0));

    testing.expect(mat.eql(transposed.transpose()));
}

test "Multiplying by a scalar" {
    const Mat1x4 = Matrix(f32, 1, 4);
    const mat = Mat1x4.fromValues([_]f32{
        7.5, 21, 0, 5,
    }, .RowMajor);

    var result = mat.multiplyScalar(2);
    const correct = Mat1x4.fromValues([_]f32{
        15, 42, 0, 10,
    }, .RowMajor);

    testing.expect(result.eql(correct));
}

test "Multiplication of two matrices" {
    const Mat1x3 = Matrix(f32, 1, 3);
    const Mat3x4 = Matrix(f32, 3, 4);

    const a = Mat1x3.fromValues([_]f32{
        3, 4, 2,
    }, .RowMajor);

    const b = Mat3x4.fromValues([_]f32{
        13, 9, 7, 15,
        8,  7, 4, 6,
        6,  4, 0, 3,
    }, .RowMajor);

    const multiplied = a.multiply(b);
    const correct = @TypeOf(multiplied).fromValues([_]f32{
        83, 63, 37, 75,
    }, .RowMajor);

    testing.expectEqual(@as(usize, 4), @TypeOf(multiplied).columns);
    testing.expectEqual(@as(usize, 1), @TypeOf(multiplied).rows);
    testing.expect(multiplied.eql(correct));
}

test "Identity Matrix" {
    const Mat4 = Matrix(u8, 4, 4);
    const Mat4x3 = Matrix(u8, 4, 3);

    testing.expect(@hasDecl(Mat4, "identity"));
    testing.expect(!@hasDecl(Mat4x3, "identity"));

    const id = Mat4.fromValues([_]u8{
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1,
    }, .RowMajor);

    testing.expect(Mat4.identity.eql(id));
}

test "Multiplication by identity" {
    const Mat4 = Matrix(u8, 4, 4);

    var id = Mat4.identity;
    var a = Mat4.fromValues([_]u8{
        1,  2,  3,  4,
        5,  6,  7,  8,
        9,  10, 11, 12,
        13, 14, 15, 16,
    }, .RowMajor);

    testing.expect(a.eql(a.multiply(id)));
    testing.expect(a.eql(id.multiply(a)));
    testing.expect(a.eql(id.multiply(a).multiply(id).multiply(id)));
}


test "fromValues" {
    const MatRM = Matrix(u8, 3, 4);

    const row_major_layout = [_]u8{
        5,  9,  4, 7,
        42, 70, 3, 0,
        6,  8,  2, 1,
    };

    const mat_rm = MatRM.fromValues([_]u8{
        5,  9,  4, 7,
        42, 70, 3, 0,
        6,  8,  2, 1,
    }, .RowMajor);

    for (mat_rm.data) |v, i| {
        testing.expectEqual(row_major_layout[i], v);
    }
}

test "Determinant of a Matrix" {
    const Mat3 = Matrix(f32, 3, 3);
    const Mat4x3 = Matrix(f32, 4, 3);

    testing.expect(@hasDecl(Mat3, "determinant"));
    testing.expect(!@hasDecl(Mat4x3, "determinant"));

    const mat = Mat3.fromValues([_]f32{
        350, 30, 20,
        37,  38, 28,
        86,  71, 37,
    }, .RowMajor);

    const id_determinant = Mat3.identity.determinant();

    testing.expectEqual(@as(f32, -185350), mat.determinant());
    testing.expectEqual(@as(f32, 1), id_determinant);
}

test "Inverse of a Matrix" {
    const Mat3 = Matrix(f64, 3, 3);
    const Mat4x3 = Matrix(f32, 4, 3);

    testing.expect(@hasDecl(Mat3, "invert"));
    testing.expect(!@hasDecl(Mat4x3, "invert"));

    const mat = Mat3.fromValues([_]f64{
        45.4, 20,      -25,
        96,   -57.564, 23.5,
        -59,  -30,     5.3,
    }, .RowMajor);

    const correct = Mat3.fromValues([_]f64{
        0.0029157187487680710108, 0.0046953542495142360015,  -0.0070656332347892020334,
        -0.013818485883702378095, -0.0089997692212971780058, -0.025276900073976560239,
        -0.045759843459199085518, 0.0013269479400801101745,  -0.033052710013558439084,
    }, .RowMajor);

    const inverse = try mat.invert();

    testing.expect(inverse.approxEql(correct, 0.00000000000001));

    const inv_mat = inverse.multiply(mat);
    const mat_inv = mat.multiply(inverse);

    testing.expect(inv_mat.approxEql(Mat3.identity, 0.00000000000001));
    testing.expect(mat_inv.approxEql(Mat3.identity, 0.00000000000001));

    const singular = Mat3.fromValues([_]f64{
        45.4, 0, -25,
        96,   0, 23.5,
        -59,  0, 5.3,
    }, .RowMajor);

    testing.expectError(error.SingularMatrix, singular.invert());
}

test "getData" {
    const Mat2 = Matrix(u8, 2, 2);
    const mat = Mat2.fromValues([_]u8{
        1, 2,
        3, 4,
    }, .RowMajor);

    const row_major = mat.getData(.RowMajor);
    const column_major = mat.getData(.ColumnMajor);

    testing.expect(std.mem.eql(u8, &row_major, &([_]u8{ 1, 2, 3, 4 })));
    testing.expect(std.mem.eql(u8, &column_major, &([_]u8{ 1, 3, 2, 4 })));
}
