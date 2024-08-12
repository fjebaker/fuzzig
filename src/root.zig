const std = @import("std");
const utils = @import("utils.zig");
const structures = @import("structures.zig");

const CharacterType = utils.CharacterType;
const MatrixT = structures.MatrixT;

pub const Unicode = if (@import("options").unicode)
    @import("unicode.zig").Unicode
else
    @compileError("Not compiled with unicode support");

test "other" {
    comptime if (@import("options").unicode) {
        _ = @import("unicode.zig");
    };
}

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

pub fn ScoresType(comptime ScoreT: type) type {
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

pub fn AlgorithmType(
    comptime ElType: type,
    comptime ScoreT: type,
    comptime scores: ScoresType(ScoreT),
) type {
    return struct {
        const Matrix = MatrixT(ScoreT);
        const Self = @This();

        pub fn FunctionTable(
            comptime Ctx: type,
        ) type {
            return struct {
                isEqual: fn (Ctx, ElType, ElType) bool,
                bonus: fn (Ctx, ScoresType(ScoreT), ElType, ElType) ScoreT,
                score: fn (Ctx, ScoresType(ScoreT), ElType, ElType) ?ScoreT,
            };
        }

        // Scoring matrix
        m: Matrix,
        // Skip score matrix
        x: Matrix,
        // Skip index
        m_skip: MatrixT(bool),

        // Character positional / role bonuses
        role_bonus: []ScoreT,
        // Buffer for the current row's bonuses
        bonus_buffer: []ScoreT,
        // Buffer for the row wise first matches
        first_match_buffer: []usize,

        traceback_buffer: []usize,

        max_haystack: usize,
        max_needle: usize,

        allocator: std.mem.Allocator,

        pub fn deinit(self: *Self) void {
            self.m.deinit();
            self.x.deinit();
            self.m_skip.deinit();
            self.allocator.free(self.role_bonus);
            self.allocator.free(self.bonus_buffer);
            self.allocator.free(self.first_match_buffer);
            self.allocator.free(self.traceback_buffer);
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
            var m_skip = try MatrixT(bool).init(allocator, rows, cols);
            errdefer m_skip.deinit();

            const role_bonus = try allocator.alloc(ScoreT, cols);
            errdefer allocator.free(role_bonus);

            const bonus_buffer = try allocator.alloc(ScoreT, cols);
            errdefer allocator.free(bonus_buffer);

            const first_match_buffer = try allocator.alloc(usize, rows);
            errdefer allocator.free(first_match_buffer);

            const traceback_buffer = try allocator.alloc(usize, cols);
            errdefer allocator.free(traceback_buffer);

            return .{
                .m = m,
                .x = x,
                .m_skip = m_skip,

                .role_bonus = role_bonus,
                .bonus_buffer = bonus_buffer,
                .first_match_buffer = first_match_buffer,
                .traceback_buffer = traceback_buffer,

                .max_haystack = max_haystack,
                .max_needle = max_needle,

                .allocator = allocator,
            };
        }

        /// Resize pre-allocated buffers to fit a new maximum haystack and
        /// needle size
        pub fn resize(self: *Self, max_haystack: usize, max_needle: usize) !void {
            const new_rows = max_needle + 1;
            const new_cols = max_haystack + 1;

            self.first_match_buffer = try self.allocator.realloc(
                self.first_match_buffer,
                new_rows,
            );

            self.role_bonus = try self.allocator.realloc(
                self.role_bonus,
                new_cols,
            );

            self.bonus_buffer = try self.allocator.realloc(
                self.bonus_buffer,
                new_cols,
            );

            self.traceback_buffer = try self.allocator.realloc(
                self.traceback_buffer,
                new_cols,
            );

            try self.m.resizeAlloc(new_rows, new_cols);
            try self.x.resizeAlloc(new_rows, new_cols);
            try self.m_skip.resizeAlloc(new_rows, new_cols);

            self.max_haystack = max_haystack;
            self.max_needle = max_needle;
        }

        /// Ensure that there is enought memory allocated, and allocate memory only if needed.
        /// Will never shrink memory.
        pub fn ensureSize(self: *Self, max_haystack: usize, max_needle: usize) !void {
            const new_rows = max_needle + 1;
            const new_cols = max_haystack + 1;

            if (new_rows < self.max_needle and new_cols < self.max_haystack) {
                return;
            }

            if (new_rows > self.max_needle + 1) {
                self.first_match_buffer = try self.allocator.realloc(
                    self.first_match_buffer,
                    new_rows,
                );
                self.max_needle = max_needle;
            }
            if (new_cols > self.max_haystack + 1) {
                self.role_bonus = try self.allocator.realloc(
                    self.role_bonus,
                    new_cols,
                );

                self.bonus_buffer = try self.allocator.realloc(
                    self.bonus_buffer,
                    new_cols,
                );

                self.traceback_buffer = try self.allocator.realloc(
                    self.traceback_buffer,
                    new_cols,
                );
                self.max_haystack = max_haystack;
            }
            if (new_cols * new_rows < (self.max_haystack + 1) * (self.max_needle + 1)) {
                try self.m.resizeAlloc(new_rows, new_cols);
                try self.x.resizeAlloc(new_rows, new_cols);
                try self.m_skip.resizeAlloc(new_rows, new_cols);
            }
            return;
        }

        /// Compute matching score
        pub fn score(
            self: *Self,
            ctx: anytype,
            comptime funcTable: FunctionTable(@TypeOf(ctx)),
            haystack: []const ElType,
            needle: []const ElType,
        ) ?ScoreT {
            const info = self.scoreImpl(ctx, funcTable, haystack, needle) orelse
                return null;
            return info.score;
        }

        /// Return the maximum (`haystack.len <= MAXIMUM`) haystack length that
        /// the algorithm has allocated memory for. To increase, use `resize`
        pub fn maximumHaystackLen(self: *const Self) usize {
            return self.max_haystack;
        }

        /// Return the maximum (`needle.len <= MAXIMUM`) needle length that
        /// the algorithm has allocated memory for. To increase, use `resize`
        pub fn maximumNeedleLen(self: *const Self) usize {
            return self.max_needle;
        }

        const ScoreInfo = struct {
            score: ScoreT,
            col_max: usize = 0,
            first_match_indices: []const usize = &.{},
        };

        fn scoreImpl(
            self: *Self,
            ctx: anytype,
            comptime funcTable: FunctionTable(@TypeOf(ctx)),
            haystack: []const ElType,
            needle: []const ElType,
        ) ?ScoreInfo {
            if (needle.len == 0) return .{
                .score = 0,
            };

            std.debug.assert(haystack.len <= self.maximumHaystackLen());
            std.debug.assert(needle.len < self.maximumNeedleLen());

            const rows = needle.len;
            const cols = haystack.len;

            // resize the view into memory
            self.m.resizeNoAlloc(rows + 1, cols + 1);
            self.x.resizeNoAlloc(rows + 1, cols + 1);
            self.m_skip.resizeNoAlloc(rows + 1, cols + 1);

            const first_match_indices = utils.firstMatchesGeneric(
                ElType,
                ctx,
                funcTable.isEqual,
                self.first_match_buffer,
                haystack,
                needle,
            ) orelse return null;

            self.reset(rows + 1, cols + 1, first_match_indices);
            self.determineBonuses(ctx, funcTable, haystack);

            self.populateMatrices(
                ctx,
                funcTable,
                haystack,
                needle,
                first_match_indices,
            );
            const col_max = self.findMaximalElement(
                first_match_indices,
                rows,
                cols,
            );

            const last_row_index = needle.len;
            const s = self.m.get(last_row_index, col_max);
            return .{
                .score = s,
                .col_max = col_max,
                .first_match_indices = first_match_indices,
            };
        }

        pub const Matches = struct {
            score: ?ScoreT,
            matches: []const usize = &.{},
        };

        /// Compute the score and the indices of the matched characters
        pub fn scoreMatches(
            self: *Self,
            ctx: anytype,
            comptime funcTable: FunctionTable(@TypeOf(ctx)),
            haystack: []const ElType,
            needle: []const ElType,
        ) Matches {
            const s = self.scoreImpl(ctx, funcTable, haystack, needle) orelse
                return .{ .score = null };

            const matches = self.traceback(
                needle.len,
                s.col_max,
                self.traceback_buffer,
                s.first_match_indices,
            );
            return .{ .score = s.score, .matches = matches };
        }

        fn traceback(
            self: *Self,
            row_max: usize,
            col_max: usize,
            traceback_buffer: []usize,
            first_match_indices: []const usize,
        ) []const usize {
            var row = row_max;
            var col = col_max;

            var tb_index: usize = 0;

            // the first visited element in the first row
            const beginning = first_match_indices[0];

            while (row > 0 and col > beginning) {
                const skipped = self.m_skip.get(row, col);
                if (skipped == false) {
                    traceback_buffer[tb_index] = col - 1;
                    tb_index += 1;
                    row -= 1;
                }
                col -= 1;
            }

            const buf = traceback_buffer[0..tb_index];
            std.mem.reverse(usize, buf);
            return buf;
        }

        fn determineBonuses(
            self: *Self,
            ctx: anytype,
            comptime funcTable: FunctionTable(@TypeOf(ctx)),
            haystack: []const ElType,
        ) void {
            var prev: ElType = 0;
            for (1.., haystack) |i, h| {
                self.role_bonus[i] = funcTable.bonus(ctx, scores, prev, h);
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
            cols: usize,
            first_match_indices: []const usize,
        ) void {
            self.m.fill(0);
            self.x.fill(0);
            self.m_skip.fill(true);
            @memset(self.first_match_buffer[0..rows], 0);
            @memset(self.traceback_buffer[0..cols], 0);

            for (1..rows) |i| {
                if (first_match_indices.len <= i - 1) break;
                const j = first_match_indices[i - 1];
                self.m.set(i, j, scores.default_score);
                self.x.set(i, j, scores.default_score);
            }

            @memset(self.role_bonus[0..cols], 0);
            @memset(self.bonus_buffer[0..cols], 0);
        }

        fn findMaximalElement(
            self: *const Self,
            first_match_indices: []const usize,
            last_row_index: usize,
            last_col_index: usize,
        ) usize {
            // the first visted element of the last row
            const start =
                first_match_indices[first_match_indices.len - 1];
            const last_row = self.m.getRow(last_row_index);

            // iterate over the visited element of the last row
            const index = std.mem.indexOfMax(
                ScoreT,
                last_row[start .. last_col_index + 1],
            );
            return index + start;
        }

        fn populateMatrices(
            self: *Self,
            ctx: anytype,
            comptime funcTable: FunctionTable(@TypeOf(ctx)),
            haystack: []const ElType,
            needle: []const ElType,
            first_match_indices: []const usize,
        ) void {
            for (1.., needle) |i, n| {

                // how many characters of the haystack do we skip
                const skip = first_match_indices[i - 1];
                for (skip + 1.., haystack[skip..]) |j, h| {
                    const prev_score_match = self.m.get(i - 1, j - 1);
                    const prev_score_skip = self.x.get(i - 1, j - 1);

                    // start by updating the M matrix

                    // compute score
                    if (funcTable.score(ctx, scores, h, n)) |current| {
                        const prev_bonus = self.bonus_buffer[j - 1];

                        // role bonus for current character
                        const role_bonus = self.role_bonus[j];

                        // TODO: this could be iffy? if too many penalties arrise
                        // then there is a chance that this is negative
                        // but that the previous index was still matched
                        // const prev_matched = prev_score_match > 0;
                        const prev_matched = !self.m_skip.get(i - 1, j - 1);
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

                        self.m_skip.set(i, j, false);

                        if (score_match >= score_skip) {
                            self.m.set(i, j, score_match);
                        } else {
                            self.m.set(i, j, score_skip);
                        }
                    } else {
                        self.m.set(i, j, scores.default_score);
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
                    } else {
                        self.x.set(i, j, x_extend);
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

            const max_col = haystack.len + 1;

            try writer.writeByteNTimes(' ', 4);
            for (self.m.getRow(0)[0..max_col], self.m_skip.getRow(0)[0..max_col]) |el, skip| {
                const width = utils.digitCount(el);
                try writer.writeByteNTimes(' ', el_width - width);
                try writer.print("{d}{c}", .{ el, if (skip) @as(u8, '-') else @as(u8, 'm') });
            }

            try writer.writeByte('\n');
            for (needle, 1..) |n, row| {
                try writer.print("'{c}'", .{n});
                try writer.writeByteNTimes(' ', 1);

                for (self.m.getRow(row)[0..max_col], self.m_skip.getRow(row)[0..max_col]) |el, skip| {
                    const width = utils.digitCount(el);
                    try writer.writeByteNTimes(' ', el_width - width);
                    try writer.print("{d}{c}", .{ el, if (skip) @as(u8, '-') else @as(u8, 'm') });
                }

                try writer.writeByte('\n');
            }
        }
    };
}

pub const Ascii = struct {
    pub const Algorithm = AlgorithmType(u8, i32, .{});
    pub const Scores = ScoresType(i32);

    const FunctionTable: Algorithm.FunctionTable(*Ascii) = .{
        .score = scoreFunc,
        .bonus = bonusFunc,
        .isEqual = eqlFunc,
    };

    fn eqlFunc(self: *Ascii, h: u8, n: u8) bool {
        if (n == ' ' and self.opts.wildcard_spaces) {
            return switch (h) {
                'a'...'z', 'A'...'Z', '0'...'9' => false,
                else => true,
            };
        } else if (!self.opts.case_sensitive) {
            return std.ascii.toLower(h) == std.ascii.toLower(n);
        } else {
            return h == n;
        }
    }

    fn scoreFunc(
        a: *Ascii,
        scores: Scores,
        h: u8,
        n: u8,
    ) ?i32 {
        if (!a.eqlFunc(h, n)) return null;

        if (a.opts.case_penalize and (h != n)) {
            return scores.score_match + a.opts.penalty_case_mistmatch;
        }
        return scores.score_match;
    }

    fn bonusFunc(
        _: *Ascii,
        scores: Scores,
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

    pub const Options = struct {
        case_sensitive: bool = true,
        case_penalize: bool = false,
        // treat spaces as wildcards for any kind of boundary
        // i.e. match with any `[^a-z,A-Z,0-9]`
        wildcard_spaces: bool = false,

        penalty_case_mistmatch: i32 = -2,
    };

    alg: Algorithm,
    opts: Options,

    // public interface

    pub fn init(
        allocator: std.mem.Allocator,
        max_haystack: usize,
        max_needle: usize,
        opts: Options,
    ) !Ascii {
        const alg = try Algorithm.init(allocator, max_haystack, max_needle);
        return Ascii{ .alg = alg, .opts = opts };
    }

    pub fn deinit(self: *Ascii) void {
        self.alg.deinit();
    }

    /// Compute matching score
    pub fn score(
        self: *Ascii,
        haystack: []const u8,
        needle: []const u8,
    ) ?i32 {
        return self.alg.score(self, FunctionTable, haystack, needle);
    }

    /// Compute the score and the indices of the matched characters
    pub fn scoreMatches(
        self: *Ascii,
        haystack: []const u8,
        needle: []const u8,
    ) Algorithm.Matches {
        return self.alg.scoreMatches(self, FunctionTable, haystack, needle);
    }

    /// Resize pre-allocated buffers to fit a new maximum haystack and
    /// needle size
    pub fn resize(self: *Ascii, max_haystack: usize, max_needle: usize) !void {
        try self.alg.resize(max_haystack, max_needle);
    }
};

fn doTestScore(alg: *Ascii, haystack: []const u8, needle: []const u8, comptime score: i32) !void {
    const s = alg.score(haystack, needle);
    try std.testing.expectEqual(score, s.?);
}

test "algorithm test" {
    const o = Ascii.Scores{};

    var alg = try Ascii.init(
        std.testing.allocator,
        128,
        32,
        .{},
    );
    defer alg.deinit();

    try doTestScore(&alg, "ab", "ab", o.score_match * 2 +
        (o.bonus_head * o.bonus_first_character_multiplier) +
        o.bonus_consecutive);

    try doTestScore(&alg, "xab", "ab", o.score_match * 2 +
        o.bonus_consecutive);

    try doTestScore(&alg, "xabbaababab", "aba", o.score_match * 3 +
        o.bonus_consecutive * 2);

    try doTestScore(&alg, "aoo_boo", "ab", o.score_match * 2 +
        (o.bonus_head * o.bonus_first_character_multiplier) +
        o.score_gap_start + o.score_gap_extension * 3 +
        o.bonus_break);

    try doTestScore(&alg, "acb", "ab", o.score_match * 2 +
        (o.bonus_head * o.bonus_first_character_multiplier) +
        o.score_gap_start + o.score_gap_extension);

    try doTestScore(&alg, "hello world", "orld", o.score_match * 4 +
        o.bonus_consecutive * 3);

    try doTestScore(&alg, "hello world", "wrld", o.score_match * 4 +
        o.bonus_consecutive * 2 +
        o.bonus_head +
        o.score_gap_start + o.score_gap_extension * 1);

    try doTestScore(&alg, "hello woorld", "wrld", o.score_match * 4 +
        o.bonus_consecutive * 2 +
        o.bonus_head +
        o.score_gap_start + o.score_gap_extension * 2);

    try doTestScore(&alg, "aaaaaaaaaaaaaaab", "b", o.score_match);

    try doTestScore(&alg, "ac" ++ "a" ** 13 ++ "b", "cb", o.score_match * 2 +
        o.score_gap_start + o.score_gap_extension * 13);

    try doTestScore(
        &alg,
        "ac" ++ "a" ** 13 ++ "b" ++ "a" ** 14 ++ "d",
        "cbd",
        o.score_match * 3 + o.score_gap_start * 2 + o.score_gap_extension * 26,
    );

    try doTestScore(&alg, "focarbaz1", "obz", o.score_match * 3 +
        o.score_gap_start * 2 +
        o.score_gap_extension * 4);

    try doTestScore(&alg, "fooBarbaz1", "obz", o.score_match * 3 +
        o.score_gap_start * 2 +
        o.score_gap_extension * 4);

    try doTestScore(&alg, "fooBarbaz1", "fobz", o.score_match * 4 +
        (o.bonus_head * o.bonus_first_character_multiplier) +
        o.bonus_consecutive +
        o.score_gap_start * 2 +
        o.score_gap_extension * 5);

    try doTestScore(&alg, "foBarbaz1", "fobz", o.score_match * 4 +
        (o.bonus_head * o.bonus_first_character_multiplier) +
        o.bonus_consecutive +
        o.score_gap_start * 2 +
        o.score_gap_extension * 4);

    try doTestScore(&alg, "fooBarbaz1", "obz", o.score_match * 3 +
        o.score_gap_start * 2 +
        o.score_gap_extension * 4);

    try doTestScore(&alg, "xhell o", "helo", o.score_match * 4 +
        o.bonus_break +
        o.bonus_consecutive * 3 +
        o.score_gap_start * 1 +
        o.score_gap_extension * 3);
}

test "case sensitivity" {
    const o = Ascii.Scores{};

    var alg1 = try Ascii.init(
        std.testing.allocator,
        128,
        32,
        .{ .case_sensitive = false },
    );
    defer alg1.deinit();

    try doTestScore(
        &alg1,
        "xab",
        "ab",
        o.score_match * 2 +
            o.bonus_consecutive,
    );

    try doTestScore(
        &alg1,
        "xaB",
        "ab",
        o.score_match * 2 +
            o.bonus_consecutive,
    );

    var alg2 = try Ascii.init(
        std.testing.allocator,
        128,
        32,
        .{
            .case_sensitive = false,
            .case_penalize = true,
        },
    );
    defer alg2.deinit();

    const A: Ascii.Options = .{};
    try doTestScore(
        &alg2,
        "xaB",
        "ab",
        o.score_match * 2 +
            o.bonus_consecutive + A.penalty_case_mistmatch,
    );
}

test "wildcard space" {
    const o = Ascii.Scores{};
    var alg = try Ascii.init(
        std.testing.allocator,
        128,
        32,
        .{ .wildcard_spaces = true },
    );
    defer alg.deinit();
    try doTestScore(
        &alg,
        "ztechno-beyond things",
        "tech",
        o.score_match * 4 +
            o.bonus_consecutive * 3,
    );

    try doTestScore(
        &alg,
        "z-abc",
        " ",
        o.score_match * 1,
    );

    try doTestScore(
        &alg,
        "ztechno-beyond things",
        "tech ",
        o.score_match * 5 +
            o.bonus_consecutive * 3 +
            o.score_gap_start + o.score_gap_extension * 2,
    );
}

fn doTestTraceback(alg: *Ascii, haystack: []const u8, needle: []const u8, comptime matches: []const usize) !void {
    const s = alg.scoreMatches(haystack, needle);

    // const stderr = std.io.getStdErr().writer();
    // try alg.debugPrint(stderr, haystack, needle);
    // std.debug.print("SCORE : {d}\n", .{s.score orelse -1});

    try std.testing.expectEqualSlices(usize, matches, s.matches);
}

test "traceback" {
    var alg = try Ascii.init(
        std.testing.allocator,
        64,
        32,
        .{},
    );
    defer alg.deinit();

    try doTestTraceback(&alg, "ab", "ab", &.{ 0, 1 });
    try doTestTraceback(&alg, "a_b", "ab", &.{ 0, 2 });
    try doTestTraceback(&alg, "abcdefg", "abcefg", &.{ 0, 1, 2, 4, 5, 6 });
    try doTestTraceback(&alg, "A" ++ "a" ** 20 ++ "B", "AB", &.{ 0, 21 });
    try doTestTraceback(&alg, "./src/main.zig", "main", &.{ 6, 7, 8, 9 });
}
