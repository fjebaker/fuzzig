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

        /// Resize the matrix without reallocating. Asserts that there is enough space.
        pub fn resizeNoAlloc(m: *Self, new_rows: usize, new_cols: usize) void {
            std.debug.assert(m.matrix.len >= new_rows * new_cols);

            m.rows = new_rows;
            m.cols = new_cols;
        }

        /// Resize the matrix to a bigger or smaller one. Note that the data will be invalided.
        pub fn resizeAlloc(m: *Self, new_rows: usize, new_cols: usize) !void {
            m.matrix = try m.allocator.realloc(m.matrix, new_cols * new_rows);
            m.resizeNoAlloc(new_rows, new_cols);
        }

        /// Set the currently active region of the matrix to a specific value.
        pub fn fill(m: Self, value: ElementType) void {
            const end = m.rows * m.cols;
            @memset(m.matrix[0..end], value);
        }
    };
}
