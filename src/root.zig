const std = @import("std");
const utils = @import("utils.zig");
const structures = @import("structures.zig");

const CharacterType = utils.CharacterType;
const MatrixT = structures.MatrixT;

// Matrix Filling: for two sequences a1a2a3..., b1b2b3... we have a reward
// function
//
//     s(a_i, b_j) = { +SCORE  if a_i == b_j     # match
//                   { -SCORE  if a_i != b_j     # mitmatch
//
//     s(-, b_j) = -PENALTY if a_i == -          # insertion
//     s(a_i, -) = -PENALTY if b_j == -          # deletion
//
// Here `-` is a placeholder to represent an "unknown" or gap in the sequence.
//
//                  { M[i-1,j-1] + s(a_i, b_i)   # up and left
//     M[i,j] = max { M[i,j-1] + s(a_i, -)       # left
//                  { M[i-1,j] + s(-, b_i)       # up
//                  { 0                          # no direction
//
// This matrix is the score of best alignment between i characters of x and j
// characters of y.
//
// For fuzzy finding want to include a gap penalty for successive gaps, i.e. so
// a run of `n` gaps incurrs less of a penalty than starting `n` individual 1
// character gaps
//
//     gaps(n) < penalty * n
//
// Gap information is contained in two new matrices tracking row and column
// gaps respectively:
//
//     X[i,j] = max { M[i,j-l] - gaps(l)     l in 1:j
//                  { Y[i,j-l] - gaps(l)     l in 1:j
//
//     Y[i,j] = max { M[i-l,j] - gaps(l)     l in 1:j
//                  { X[i-l,j] - gaps(l)     l in 1:j
//
// Cost of building these matrices is O(n^3). Instead, assume an _affine gap
// penalty_, where we differentiate the penalty for `gap_start` and
// `gap_ext` for starting a new gap or extending the existing gap.
// Requires the same three matrices but no longer search over all gap lengths,
// only need to know if we're starting a new gap
//
//                  { gap_start + gap_ext + M[i,j-1]
//     X[i,j] = max { gap_ext + X[i,j-1]
//                  { gap_start + gap_ext + Y[i,j-1]
//
// The algorithm is now O(mn).

pub fn Scores(comptime ScoreT: type) type {
    return struct {
        score_match: ScoreT = 16,
        score_gap_start: ScoreT = -3,
        score_gap_extension: ScoreT = -1,

        default_score: ScoreT = 0,

        // bonuses
        bonus_consecutive: ScoreT = 4,
        bonus_first_character_multiplier: ScoreT = 2,

        bonus_break: ScoreT = 5,
        bonus_camel: ScoreT = 2,
        bonus_head: ScoreT = 8,
        bonus_tail: ScoreT = 0,
    };
}

pub fn Algorithm(
    comptime ElType: type,
    comptime ScoreT: type,
    comptime scores: Scores(ScoreT),
    comptime Impl: type,
) type {
    return struct {
        const Matrix = MatrixT(ScoreT);
        const Self = @This();

        // Scoring matrix
        m: Matrix,
        // Skip score matrix
        x: Matrix,
        // Skip index
        skip: MatrixT(bool),

        // Character positional / role bonuses
        role_bonus: []ScoreT,
        // Buffer for the current row's bonuses
        bonus_buffer: []ScoreT,
        // Buffer for the row wise first matches
        first_match_buffer: []usize,

        allocator: std.mem.Allocator,

        impl: Impl,

        pub fn deinit(self: *Self) void {
            self.m.deinit();
            self.x.deinit();
            self.skip.deinit();
            self.allocator.free(self.role_bonus);
            self.allocator.free(self.bonus_buffer);
            self.allocator.free(self.first_match_buffer);
            self.* = undefined;
        }

        pub fn init(
            allocator: std.mem.Allocator,
            max_haystack: usize,
            max_needle: usize,
            impl: Impl,
        ) !Self {
            const rows = max_needle + 1;
            const cols = max_haystack + 1;

            var m = try Matrix.init(allocator, rows, cols);
            errdefer m.deinit();
            var x = try Matrix.init(allocator, rows, cols);
            errdefer x.deinit();
            var skip = try MatrixT(bool).init(allocator, rows, cols);
            errdefer skip.deinit();

            const role_bonus = try allocator.alloc(ScoreT, cols);
            errdefer allocator.free(role_bonus);

            const bonus_buffer = try allocator.alloc(ScoreT, cols);
            errdefer allocator.free(bonus_buffer);

            const first_match_buffer = try allocator.alloc(usize, rows);
            errdefer allocator.free(first_match_buffer);

            return .{
                .m = m,
                .x = x,
                .skip = skip,

                .role_bonus = role_bonus,
                .bonus_buffer = bonus_buffer,
                .first_match_buffer = first_match_buffer,
                .allocator = allocator,
                .impl = impl,
            };
        }

        /// Compute matching score
        pub fn score(
            self: *Self,
            haystack: []const ElType,
            needle: []const ElType,
        ) ?ScoreT {
            if (needle.len == 0) return 0;

            const rows = needle.len;
            const cols = haystack.len;

            std.debug.assert(rows <= self.m.rows);
            std.debug.assert(cols <= self.m.cols);

            const first_match_indices = utils.firstMatchesGeneric(
                ElType,
                &self.impl,
                Impl.eqlFunc,
                self.first_match_buffer,
                haystack,
                needle,
            ) orelse return null;

            self.reset(rows, first_match_indices);
            self.determineBonuses(haystack);

            try self.populateMatrices(haystack, needle, first_match_indices);
            const col_max = self.findMaximalElement(first_match_indices, needle.len);

            // self.traceback(col_max, first_match_indices);

            const last_row_index = needle.len;
            return self.m.get(last_row_index, col_max);
        }

        fn determineBonuses(self: *Self, haystack: []const ElType) void {
            var prev: u8 = 0;
            for (1.., haystack) |i, h| {
                self.role_bonus[i] = Impl.bonusFunc(&self.impl, scores, prev, h);
                prev = h;
            }

            if (self.role_bonus.len > 1) {
                self.role_bonus[1] = self.role_bonus[1] *
                    scores.bonus_first_character_multiplier;
            }
        }

        fn reset(
            self: *Self,
            rows: usize,
            first_match_indices: []const usize,
        ) void {
            // TODO: remove, this is for debugging
            @memset(self.m.matrix, 0);
            @memset(self.x.matrix, 0);
            @memset(self.skip.matrix, true);

            // set the first row and column to zero
            @memset(self.m.getRow(0), scores.default_score);
            @memset(self.x.getRow(0), scores.default_score);

            for (1..rows) |i| {
                if (first_match_indices.len <= i - 1) break;
                const j = first_match_indices[i - 1];
                self.m.set(i, j, scores.default_score);
                self.x.set(i, j, scores.default_score);
            }

            @memset(self.bonus_buffer, 0);
            @memset(self.role_bonus, 0);
        }

        fn findMaximalElement(
            self: *const Self,
            first_match_indices: []const usize,
            last_row_index: usize,
        ) usize {
            // the first visted element of the last row
            const start =
                first_match_indices[first_match_indices.len - 1];
            const last_row = self.m.getRow(last_row_index);

            // iterate over the visited element of the last row
            return std.mem.indexOfMax(ScoreT, last_row[start..]) + start;
        }

        fn populateMatrices(
            self: *Self,
            haystack: []const ElType,
            needle: []const ElType,
            first_match_indices: []const usize,
        ) !void {
            for (1.., needle) |i, n| {

                // how many characters of the haystack do we skip
                const skip = first_match_indices[i - 1];
                for (skip + 1.., haystack[skip..]) |j, h| {
                    const prev_score_match = self.m.get(i - 1, j - 1);
                    const prev_score_skip = self.x.get(i - 1, j - 1);

                    // start by updating the M matrix

                    // compute score
                    if (Impl.scoreFunc(&self.impl, scores, h, n)) |current| {
                        const prev_bonus = self.bonus_buffer[j - 1];

                        // role bonus for current character
                        const role_bonus = self.role_bonus[j];

                        const prev_matched = !self.skip.get(i - 1, j - 1);
                        const consecutive_bonus = if (prev_matched)
                            scores.bonus_consecutive
                        else
                            0;

                        const current_bonus = @max(
                            prev_bonus,
                            @max(role_bonus, consecutive_bonus),
                        );

                        self.bonus_buffer[j - 1] = current_bonus;

                        const score_match = prev_score_match + current_bonus + current;
                        const score_skip = prev_score_skip + role_bonus + current;

                        if (score_match >= score_skip) {
                            self.m.set(i, j, score_match);
                            self.skip.set(i, j - 1, false);
                        } else {
                            self.m.set(i, j, score_skip);
                            self.skip.set(i, j - 1, true);
                        }
                    } else {
                        self.m.set(i, j, scores.default_score);
                        self.skip.set(i, j - 1, true);
                        self.bonus_buffer[j] = 0;
                    }

                    // then update the X matrix

                    // cost of starting a new gap
                    const x_start =
                        scores.score_gap_start +
                        scores.score_gap_extension +
                        self.m.get(i, j - 1); // M[i,j-1]

                    // cost of extending
                    const x_extend =
                        scores.score_gap_extension +
                        self.x.get(i, j - 1); // X[i,j-1]

                    if (x_start >= x_extend) {
                        self.x.set(i, j, x_start);
                        self.skip.set(i, j - 1, false);
                    } else {
                        self.x.set(i, j, x_extend);
                        self.skip.set(i, j - 1, true);
                    }
                }
            }
        }

        fn debugPrint(
            self: *const Self,
            writer: anytype,
            haystack: []const ElType,
            needle: []const ElType,
        ) !void {
            const el_width = bonus: {
                var max_digits: usize = 1;
                for (self.m.matrix) |i| {
                    max_digits = @max(utils.digitCount(i), max_digits);
                }
                break :bonus max_digits + 3;
            };

            // padding
            try writer.writeByteNTimes(' ', el_width + 3);

            // upper right corner
            try writer.writeByteNTimes(' ', el_width);
            for (haystack) |h| {
                try writer.print("'{c}'", .{h});
                try writer.writeByteNTimes(' ', el_width - 2);
            }

            try writer.writeByteNTimes('\n', 2);

            try writer.writeByteNTimes(' ', 4);
            for (self.m.getRow(0), self.skip.getRow(0)) |el, skip| {
                const width = utils.digitCount(el);
                try writer.writeByteNTimes(' ', el_width - width);
                try writer.print("{d}{c}", .{ el, if (skip) @as(u8, '-') else @as(u8, 'm') });
            }

            try writer.writeByte('\n');
            for (needle, 1..) |n, row| {
                try writer.print("'{c}'", .{n});
                try writer.writeByteNTimes(' ', 1);

                for (self.m.getRow(row), self.skip.getRow(row)) |el, skip| {
                    const width = utils.digitCount(el);
                    try writer.writeByteNTimes(' ', el_width - width);
                    try writer.print("{d}{c}", .{ el, if (skip) @as(u8, '-') else @as(u8, 'm') });
                }

                try writer.writeByte('\n');
            }
        }
    };
}

pub const AsciiOptions = struct {
    const AsciiScores = Scores(i32);

    case_sensitive: bool = true,
    case_penalize: bool = false,
    // treat spaces as wildcards for any kind of boundary
    // i.e. match with any `[^a-z,A-Z,0-9]`
    wildcard_spaces: bool = false,

    penalty_case_mistmatch: i32 = -2,

    fn eqlFunc(a: *const AsciiOptions, h: u8, n: u8) bool {
        if (n == ' ' and a.wildcard_spaces) {
            return switch (h) {
                'a'...'z', 'A'...'Z', '0'...'9' => false,
                else => true,
            };
        } else if (!a.case_sensitive) {
            return std.ascii.toLower(h) == std.ascii.toLower(n);
        } else {
            return h == n;
        }
    }

    fn scoreFunc(
        a: *const AsciiOptions,
        comptime scores: AsciiScores,
        h: u8,
        n: u8,
    ) ?i32 {
        if (!a.eqlFunc(h, n)) return null;

        if (a.case_penalize and (h != n)) {
            return scores.score_match + a.penalty_case_mistmatch;
        }
        return scores.score_match;
    }

    fn bonusFunc(
        _: *const AsciiOptions,
        comptime scores: AsciiScores,
        h: u8,
        n: u8,
    ) i32 {
        const p = CharacterType.fromAscii(h);
        const c = CharacterType.fromAscii(n);

        return switch (p.roleNextTo(c)) {
            .Head => scores.bonus_head,
            .Camel => scores.bonus_camel,
            .Break => scores.bonus_break,
            .Tail => scores.bonus_tail,
        };
    }
};

/// Default ASCII Fuzzy Finder
pub const Ascii = Algorithm(u8, i32, .{}, AsciiOptions);

fn doTestScore(opts: AsciiOptions, haystack: []const u8, needle: []const u8, comptime score: i32) !void {
    var alg = try Ascii.init(
        std.testing.allocator,
        haystack.len,
        needle.len,
        opts,
    );
    defer alg.deinit();

    const s = alg.score(haystack, needle);

    // const stderr = std.io.getStdErr().writer();
    // try alg.debugPrint(stderr, haystack, needle);
    // std.debug.print("SCORE : {d}\n", .{s orelse -1});

    try std.testing.expectEqual(score, s.?);
}

test "algorithm test" {
    const o = AsciiOptions.AsciiScores{};
    try doTestScore(.{}, "ab", "ab", o.score_match * 2 +
        (o.bonus_head * o.bonus_first_character_multiplier) +
        o.bonus_consecutive);

    try doTestScore(.{}, "xab", "ab", o.score_match * 2 +
        o.bonus_consecutive);

    try doTestScore(.{}, "xabbaababab", "aba", o.score_match * 3 +
        o.bonus_consecutive * 2);

    try doTestScore(.{}, "aoo_boo", "ab", o.score_match * 2 +
        (o.bonus_head * o.bonus_first_character_multiplier) +
        o.score_gap_start + o.score_gap_extension * 3 +
        o.bonus_break);

    try doTestScore(.{}, "acb", "ab", o.score_match * 2 +
        (o.bonus_head * o.bonus_first_character_multiplier) +
        o.score_gap_start + o.score_gap_extension);

    try doTestScore(.{}, "hello world", "orld", o.score_match * 4 +
        o.bonus_consecutive * 3);

    try doTestScore(.{}, "hello world", "wrld", o.score_match * 4 +
        o.bonus_consecutive * 2 +
        o.bonus_head +
        o.score_gap_start + o.score_gap_extension * 1);

    try doTestScore(.{}, "hello woorld", "wrld", o.score_match * 4 +
        o.bonus_consecutive * 2 +
        o.bonus_head +
        o.score_gap_start + o.score_gap_extension * 2);

    try doTestScore(.{}, "aaaaaaaaaaaaaaab", "b", o.score_match);
    try doTestScore(.{}, "acaaaaaaaaaaaaab", "cb", o.score_match * 2 +
        o.score_gap_start + o.score_gap_extension * 13);

    try doTestScore(.{}, "focarbaz1", "obz", o.score_match * 3 +
        o.score_gap_start * 2 +
        o.score_gap_extension * 4);

    try doTestScore(.{}, "fooBarbaz1", "obz", o.score_match * 3 +
        o.score_gap_start * 2 +
        o.score_gap_extension * 4);

    try doTestScore(.{}, "fooBarbaz1", "fobz", o.score_match * 4 +
        (o.bonus_head * o.bonus_first_character_multiplier) +
        o.bonus_consecutive +
        o.score_gap_start * 2 +
        o.score_gap_extension * 5);

    try doTestScore(.{}, "foBarbaz1", "fobz", o.score_match * 4 +
        (o.bonus_head * o.bonus_first_character_multiplier) +
        o.bonus_consecutive +
        o.score_gap_start * 2 +
        o.score_gap_extension * 4);

    try doTestScore(.{}, "fooBarbaz1", "obz", o.score_match * 3 +
        o.score_gap_start * 2 +
        o.score_gap_extension * 4);

    try doTestScore(.{}, "xhell o", "helo", o.score_match * 4 +
        o.bonus_break +
        o.bonus_consecutive * 3 +
        o.score_gap_start * 1 +
        o.score_gap_extension * 3);
}

test "case sensitivity" {
    const o = AsciiOptions.AsciiScores{};

    try doTestScore(
        .{ .case_sensitive = true },
        "xab",
        "ab",
        o.score_match * 2 +
            o.bonus_consecutive,
    );

    try doTestScore(
        .{ .case_sensitive = false },
        "xaB",
        "ab",
        o.score_match * 2 +
            o.bonus_consecutive,
    );

    const opts: AsciiOptions = .{
        .case_sensitive = false,
        .case_penalize = true,
    };
    try doTestScore(
        opts,
        "xaB",
        "ab",
        o.score_match * 2 +
            o.bonus_consecutive + opts.penalty_case_mistmatch,
    );
}
