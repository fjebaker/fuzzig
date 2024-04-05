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
// We no longer require the Y matrix, and the algorithm is now O(mn).

pub fn Options(comptime ScoreT: type) type {
    return struct {
        score_match: ScoreT = 16,
        score_gap_start: ScoreT = -3,
        score_gap_extension: ScoreT = -1,

        default_score: ScoreT = 0,

        // bonuses
        bonus_first_character_multiplier: ScoreT = 2,
        bonus_head: ScoreT = 8,
        bonus_camel: ScoreT = 2,
        bonus_break: ScoreT = 5,
        bonus_tail: ScoreT = 0,
        bonus_consecutive: ScoreT = 4,
    };
}

pub fn Algorithm(
    comptime ElType: type,
    comptime ScoreT: type,
    comptime opts: Options(ScoreT),
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
            std.debug.print("rows {d}/{d} cols {d}/{d}\n", .{
                rows,
                self.m.rows,
                cols,
                self.m.cols,
            });

            std.debug.assert(rows <= self.m.rows);
            std.debug.assert(cols <= self.m.cols);

            const first_match_indices = utils.firstMatches(
                ElType,
                self.first_match_buffer,
                haystack,
                needle,
            ) orelse return null;

            self.reset(rows, first_match_indices);
            self.determineBonuses(haystack);

            std.debug.print("BONUS: {any}\n", .{self.role_bonus});

            try self.populateMatrices(haystack, needle, first_match_indices);
            const col_max = self.findMaximalElement(first_match_indices, needle.len);

            // self.traceback(col_max, first_match_indices);

            const last_row_index = needle.len;
            return self.m.get(last_row_index, col_max);
        }

        fn determineBonuses(self: *Self, haystack: []const ElType) void {
            var prev: u8 = 0;
            for (1.., haystack) |i, h| {
                self.role_bonus[i] = bonusFor(prev, h);
                prev = h;
            }

            if (self.role_bonus.len > 1) {
                self.role_bonus[1] = self.role_bonus[1] *
                    opts.bonus_first_character_multiplier;
            }
        }

        fn scoreFunc(n: ElType, h: ElType) ?ScoreT {
            if (n != h) return null;
            return opts.score_match;
        }

        fn bonusFor(prev: ElType, curr: ElType) ScoreT {
            const p = CharacterType.fromAscii(prev);
            const c = CharacterType.fromAscii(curr);

            return switch (p.roleNextTo(c)) {
                .Head => opts.bonus_head,
                .Camel => opts.bonus_camel,
                .Break => opts.bonus_break,
                .Tail => opts.bonus_tail,
            };
        }

        fn reset(
            self: *Self,
            rows: usize,
            first_match_indices: []const usize,
        ) void {
            // TODO: remove, this is for debugging
            @memset(self.m.matrix, 0);
            @memset(self.x.matrix, 0);

            // set the first row and column to zero
            @memset(self.m.getRow(0), opts.default_score);
            @memset(self.x.getRow(0), opts.default_score);

            for (1..rows) |i| {
                if (first_match_indices.len <= i - 1) break;
                const j = first_match_indices[i - 1];
                self.m.set(i, j, opts.default_score);
                self.x.set(i, j, opts.default_score);
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
                    if (scoreFunc(n, h)) |current| {
                        const prev_bonus = self.bonus_buffer[j - 1];
                        const role_bonus = self.role_bonus[j];

                        const consecutive_bonus = @max(
                            prev_bonus,
                            @max(role_bonus, opts.bonus_consecutive),
                        );

                        self.bonus_buffer[j - 1] = consecutive_bonus;

                        const score_match = prev_score_match + consecutive_bonus;
                        const score_skip = prev_score_skip + role_bonus;

                        if (score_match >= score_skip) {
                            const total = current + score_match;
                            self.m.set(i, j, total);
                            self.skip.set(i, j, false);
                        } else {
                            const total = current + score_skip;
                            self.m.set(i, j, total);
                            self.skip.set(i, j, true);
                        }
                    } else {
                        self.m.set(i, j, opts.default_score);
                        self.skip.set(i, j, true);
                        self.bonus_buffer[j] = 0;
                    }

                    // then update the X matrix

                    // cost of starting a new gap
                    const x_start =
                        opts.score_gap_start +
                        opts.score_gap_extension +
                        self.m.get(i, j - 1); // M[i,j-1]

                    // cost of extending
                    const x_extend =
                        opts.score_gap_extension +
                        self.x.get(i, j - 1); // X[i,j-1]

                    if (x_start >= x_extend) {
                        self.x.set(i, j, x_start);
                        self.skip.set(i, j, false);
                    } else {
                        self.x.set(i, j, x_extend);
                        self.skip.set(i, j, true);
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
                break :bonus max_digits + 2;
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
            for (self.m.getRow(0)) |el| {
                const width = utils.digitCount(el);
                try writer.writeByteNTimes(' ', el_width - width + 1);
                try writer.print("{d}", .{el});
            }

            try writer.writeByte('\n');
            for (needle, 1..) |n, row| {
                try writer.print("'{c}'", .{n});
                try writer.writeByteNTimes(' ', 1);

                for (self.m.getRow(row)) |el| {
                    const width = utils.digitCount(el);
                    try writer.writeByteNTimes(' ', el_width - width + 1);
                    try writer.print("{d}", .{el});
                }

                try writer.writeByte('\n');
            }
        }
    };
}

pub const Ascii = Algorithm(u8, i32, .{});

fn doTest(haystack: []const u8, needle: []const u8) !void {
    var alg = try Ascii.init(
        std.testing.allocator,
        haystack.len,
        needle.len,
    );
    defer alg.deinit();

    const s = alg.score(haystack, needle);

    const stderr = std.io.getStdErr().writer();
    try alg.debugPrint(stderr, haystack, needle);

    std.debug.print("SCORE : {d}\n", .{s orelse -1});
}

fn doTestScore(haystack: []const u8, needle: []const u8, comptime score: i32) !void {
    var alg = try Ascii.init(
        std.testing.allocator,
        haystack.len,
        needle.len,
    );
    defer alg.deinit();

    const s = alg.score(haystack, needle);

    const stderr = std.io.getStdErr().writer();
    try alg.debugPrint(stderr, haystack, needle);

    std.debug.print("SCORE : {d}\n", .{s orelse -1});
    try std.testing.expectEqual(score, s.?);
}

test "algorithm test" {
    const o = Options(i32){};
    try doTestScore("ab", "ab", o.score_match * o.bonus_first_character_multiplier +
        o.score_match + o.bonus_consecutive);
    try doTest(
        "aoo_boo",
        "ab",
    );
    try doTest(
        "acb",
        "ab",
    );
    try doTest(
        "hello world",
        "orld",
    );
    try doTest(
        "hello world",
        "wrld",
    );
    try doTestScore("foBarbaz1", "obz", o.score_match * 3 +
        o.score_gap_start * 2 +
        o.score_gap_extension * 2);
    try doTestScore("fooBarbaz1", "obz", o.score_match * 3 +
        o.score_gap_start * 2 +
        o.score_gap_extension * 2);
    try doTestScore("fooBarbaz1", "fobz", o.score_match * o.bonus_first_character_multiplier +
        o.score_match * 3 +
        o.score_gap_start * 2 +
        o.score_gap_extension * 3);
    try doTestScore("foBarbaz1", "fobz", o.score_match * o.bonus_first_character_multiplier +
        o.score_match * 3 +
        o.score_gap_start * 2 +
        o.score_gap_extension * 2);
    try doTestScore("fooBarbaz1", "obz", o.score_match * 3 +
        o.score_gap_start * 2 +
        o.score_gap_extension * 2);
    try doTest(
        "foBarbaz1",
        "obz",
    );
    // try doTest(
    //     "hello world",
    //     "hello world",
    // );
}
