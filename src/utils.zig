const std = @import("std");

pub fn digitCount(v: anytype) usize {
    const abs: u32 = @intCast(@abs(v));
    if (abs == 0) return 1;
    const width: usize = @intCast(std.math.log10_int(abs));
    if (v < 0) return width + 2;
    return width + 1;
}

pub const CharacterType = enum {
    Empty,
    Upper,
    Lower,
    Number,
    HardSeperator,
    SoftSeperator,

    /// Convert an ASCII character to a `CharacterType`
    pub fn fromAscii(c: u8) CharacterType {
        return switch (c) {
            'a'...'z' => .Lower,
            'A'...'Z' => .Upper,
            '0'...'9' => .Number,
            ' ', '\\', '/', '|', '(', ')', '[', ']', '{', '}' => .HardSeperator,
            '!', '*'...'.', ':'...'@', '^'...'`', '~' => .SoftSeperator,
            0 => .Empty,
            else => .Lower,
        };
    }

    const Role = enum {
        Head,
        Break,
        Camel,
        Tail,
    };

    /// Get the `Role` of the current character positioned next to another
    /// character
    pub fn roleNextTo(s: CharacterType, o: CharacterType) Role {
        return switch (s) {
            .Empty, .HardSeperator => .Head,
            .SoftSeperator => .Break,
            .Lower => if (o == .Upper) .Camel else .Tail,
            .Number => if (o == .Upper) .Camel else .Tail,
            else => .Tail,
        };
    }
};

/// Computes a simple equality check, recording the successive indices of haystack
/// that match successive characters in needle.
pub fn firstMatches(
    comptime T: type,
    indices: []usize,
    haystack: []const T,
    needle: []const T,
) ?[]const usize {
    if (needle.len == 0) {
        return &.{};
    }
    if (needle.len > haystack.len) {
        return null;
    }

    var index: usize = 0;
    for (0.., haystack) |i, h| {
        const n = needle[index];

        if (h == n) {
            indices[index] = i;
            index += 1;
            if (index >= needle.len) break;
        }
    } else return null;

    return indices[0..index];
}

pub fn firstMatchesAlloc(
    comptime T: type,
    allocator: std.mem.Allocator,
    haystack: []const T,
    needle: []const T,
) !?[]const usize {
    const indices = try allocator.alloc(usize, needle.len);
    errdefer allocator.free(indices);
    return firstMatches(T, indices, haystack, needle);
}

fn testFirstMatch(
    haystack: []const u8,
    needle: []const u8,
    comptime expected: []const usize,
) !void {
    const inds = try firstMatchesAlloc(u8, std.testing.allocator, haystack, needle);
    defer if (inds) |x| std.testing.allocator.free(x);

    try std.testing.expectEqualSlices(usize, expected, inds.?);
}

test "firstMatches" {
    try testFirstMatch("axbycz", "xyz", &.{ 1, 3, 5 });
    try testFirstMatch("axbycz", "abc", &.{ 0, 2, 4 });
    try testFirstMatch("", "", &.{});
}
