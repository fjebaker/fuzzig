const std = @import("std");

const GenCatData = @import("GenCatData");
const CaseData = @import("CaseData");

const UnicodeToolBox = @import("root.zig").UnicodeToolBox;

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

    pub fn fromUnicode(c: u21, unicode_toolbox: UnicodeToolBox) CharacterType {
        if (unicode_toolbox.cd.isLower(c)) {
            return .Lower;
        } else if (unicode_toolbox.cd.isUpper(c)) {
            return .Upper;
        } else if (unicode_toolbox.gcd.isNumber(c)) {
            return .Number;
        } else if (switch (c) {
            ' ', '\\', '/', '|', '(', ')', '[', ']', '{', '}' => true,
            else => false,
        }) {
            return .HardSeperator;
        } else if (unicode_toolbox.gcd.isSeparator(c)) {
            return .HardSeperator;
        } else if (unicode_toolbox.gcd.isPunctuation(c) or unicode_toolbox.gcd.isSymbol(c) or unicode_toolbox.gcd.isMark(c)) {
            return .SoftSeperator;
        } else if (unicode_toolbox.gcd.isControl(c)) {
            return .Empty;
        } else {
            return .Lower; // Maybe .Empty instead ?
        }
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

pub fn firstMatchesGeneric(
    comptime T: type,
    ctx: anytype,
    comptime eqlFunc: fn (@TypeOf(ctx), h: T, n: T) bool,
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

        if (eqlFunc(ctx, h, n)) {
            indices[index] = i;
            index += 1;
            if (index >= needle.len) break;
        }
    } else return null;

    return indices[0..index];
}

fn simpleEql(comptime T: type) fn (void, T, T, std.mem.Allocator) bool {
    return struct {
        fn f(_: void, h: T, n: T, _: std.mem.Allocator) bool {
            return h == n;
        }
    }.f;
}

/// Computes a simple equality check, recording the successive indices of haystack
/// that match successive characters in needle.
pub fn firstMatches(
    comptime T: type,
    indices: []usize,
    haystack: []const T,
    needle: []const T,
    allocator: std.mem.Allocator,
) ?[]const usize {
    return firstMatchesGeneric(T, {}, simpleEql(T), indices, haystack, needle, allocator);
}

pub fn firstMatchesAlloc(
    comptime T: type,
    allocator: std.mem.Allocator,
    haystack: []const T,
    needle: []const T,
) !?[]const usize {
    const indices = try allocator.alloc(usize, needle.len);
    errdefer allocator.free(indices);
    return firstMatches(T, indices, haystack, needle, allocator);
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
