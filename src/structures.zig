const std = @import("std");

/// Column major matrix type.
pub fn MatrixT(comptime T: type) type {
    return struct {
        pub const ElementType = T;
        const Self = @This();

        allocator: std.mem.Allocator,
        matrix: []ElementType,
        rows: usize,
        cols: usize,

        pub fn init(allocator: std.mem.Allocator, rows: usize, cols: usize) !Self {
            const mem = try allocator.alloc(ElementType, cols * rows);
            return .{
                .allocator = allocator,
                .matrix = mem,
                .rows = rows,
                .cols = cols,
            };
        }

        pub fn deinit(m: *Self) void {
            m.allocator.free(m.matrix);
            m.* = undefined;
        }

        /// Get the `row` and `col` element of the matrix.
        pub fn get(m: Self, row: usize, col: usize) ElementType {
            return m.getPtr(row, col).*;
        }

        /// Set the `row` and `col` element of the matrix.
        pub fn set(m: Self, row: usize, col: usize, t: ElementType) void {
            m.getPtr(row, col).* = t;
        }

        /// Get a pointer to `row` and `col` element of the matrix.
        pub fn getPtr(m: Self, row: usize, col: usize) *ElementType {
            return &m.matrix[m.getIndex(row, col)];
        }

        /// Get the index into the array
        pub fn getIndex(m: Self, row: usize, col: usize) usize {
            std.debug.assert(row < m.rows);
            std.debug.assert(col < m.cols);
            return row * m.cols + col;
        }

        /// Get a slice into a `row` of the matrix.
        pub fn getRow(m: Self, row: usize) []ElementType {
            std.debug.assert(row < m.rows);
            const i = row * m.cols;
            return m.matrix[i .. i + m.cols];
        }
    };
}
