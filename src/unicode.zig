const std = @import("std");
const root = @import("root.zig");
const utils = @import("utils.zig");
const structures = @import("structures.zig");

const code_point = @import("code_point");
const GenCatData = @import("GenCatData");
const CaseData = @import("CaseData");
const Normalize = @import("Normalize");

const CharacterType = utils.CharacterType;
const MatrixT = structures.MatrixT;

const AlgorithmType = root.AlgorithmType;
const ScoresType = root.ScoresType;

fn charFromUnicode(c: u21, cd: CaseData, gcd: GenCatData) CharacterType {
    if (cd.isLower(c)) {
        return .Lower;
    } else if (cd.isUpper(c)) {
        return .Upper;
    } else if (gcd.isNumber(c)) {
        return .Number;
    } else if (switch (c) {
        ' ', '\\', '/', '|', '(', ')', '[', ']', '{', '}' => true,
        else => false,
    }) {
        return .HardSeperator;
    } else if (gcd.isSeparator(c)) {
        return .HardSeperator;
    } else if (gcd.isPunctuation(c) or gcd.isSymbol(c) or gcd.isMark(c)) {
        return .SoftSeperator;
    } else if (gcd.isControl(c)) {
        return .Empty;
    } else {
        return .Lower; // Maybe .Empty instead ?
    }
}

pub const Unicode = struct {
    pub const Algorithm = AlgorithmType(u21, i32);
    pub const Scores = ScoresType(i32);

    const FunctionTable: Algorithm.FunctionTable(*Unicode) = .{
        .score = scoreFunc,
        .bonus = bonusFunc,
        .isEqual = eqlFunc,
    };

    fn eqlFunc(self: *Unicode, h: u21, n: u21) bool {
        if (self.gcd.isSeparator(n) and self.opts.wildcard_spaces) {
            if (self.gcd.isLetter(h) or
                self.gcd.isNumber(h) or
                self.gcd.isSymbol(h))
            {
                return true;
            } else {
                return false;
            }
        } else if (!self.opts.case_sensitive) {
            return self.cd.toLower(h) == self.cd.toLower(n);
        } else {
            return h == n;
        }
    }

    fn scoreFunc(
        a: *Unicode,
        scores: Scores,
        h: u21,
        n: u21,
    ) ?i32 {
        if (!a.eqlFunc(h, n)) return null;

        if (a.opts.case_penalize and (h != n)) {
            return scores.score_match + a.opts.penalty_case_mistmatch;
        }
        return scores.score_match;
    }

    fn bonusFunc(
        self: *Unicode,
        scores: Scores,
        h: u21,
        n: u21,
    ) i32 {
        const p = charFromUnicode(h, self.cd, self.gcd);
        const c = charFromUnicode(n, self.cd, self.gcd);

        return switch (p.roleNextTo(c)) {
            .Head => scores.bonus_head,
            .Camel => scores.bonus_camel,
            .Break => scores.bonus_break,
            .Tail => scores.bonus_tail,
        };
    }

    fn convertString(self: *const Unicode, string: []const u8) ![]const u21 {
        const nfc_result = try self.norm.nfc(self.alg.allocator, string);
        defer nfc_result.deinit();

        var iter = code_point.Iterator{ .bytes = nfc_result.slice };

        var converted_string = std.ArrayList(u21).init(self.alg.allocator);
        defer converted_string.deinit();

        while (iter.next()) |c| {
            try converted_string.append(c.code);
        }
        return converted_string.toOwnedSlice();
    }

    pub const Options = struct {
        case_sensitive: bool = true,
        case_penalize: bool = false,
        // treat spaces as wildcards for any kind of boundary
        // i.e. match with any `[^a-z,A-Z,0-9]`
        wildcard_spaces: bool = false,

        penalty_case_mistmatch: i32 = -2,

        char_buffer_size: usize = 8192,

        scores: Scores = .{},
    };

    alg: Algorithm,
    opts: Options,
    // unicode specific things
    gcd: GenCatData,
    norm: Normalize,
    norm_data_ptr: *Normalize.NormData,
    cd: CaseData,

    pub fn init(
        allocator: std.mem.Allocator,
        max_haystack: usize,
        max_needle: usize,
        opts: Options,
    ) !Unicode {
        var alg = try Algorithm.init(allocator, max_haystack, max_needle, opts.scores);
        errdefer alg.deinit();

        var gcd = try GenCatData.init(allocator);
        errdefer gcd.deinit();

        var norm: Normalize = undefined;
        const norm_data_ptr = try allocator.create(Normalize.NormData);
        errdefer allocator.destroy(norm_data_ptr);

        try Normalize.NormData.init(norm_data_ptr, allocator);
        errdefer norm_data_ptr.deinit();

        norm.norm_data = norm_data_ptr;

        var cd = try CaseData.init(allocator);
        errdefer cd.deinit();

        return .{
            .alg = alg,
            .opts = opts,
            .gcd = gcd,
            .norm = norm,
            .norm_data_ptr = norm_data_ptr,
            .cd = cd,
        };
    }

    pub fn deinit(self: *Unicode) void {
        self.cd.deinit();
        self.gcd.deinit();
        self.norm_data_ptr.deinit();
        self.alg.allocator.destroy(self.norm_data_ptr);
        self.alg.deinit();
    }

    /// Compute matching score. Recasts the `u8` array to `u21` to properly
    /// encode unicode characters
    pub fn score(
        self: *Unicode,
        haystack: []const u8,
        needle: []const u8,
    ) !?i32 {
        const haystack_normal = try self.convertString(haystack);
        defer self.alg.allocator.free(haystack_normal);

        const needle_normal = try self.convertString(needle);
        defer self.alg.allocator.free(needle_normal);

        return self.alg.score(
            self,
            FunctionTable,
            haystack_normal,
            needle_normal,
        );
    }

    /// Compute the score and the indices of the matched characters. Recasts
    /// the `u8` array to `u21` to properly encode unicode characters
    pub fn scoreMatches(
        self: *Unicode,
        haystack: []const u8,
        needle: []const u8,
    ) Algorithm.Matches {
        const haystack_normal = self.convertString(haystack);
        defer self.allocator.free(haystack_normal);

        const needle_normal = self.convertString(needle);
        defer self.allocator.free(needle_normal);

        return self.alg.scoreMatches(
            self,
            FunctionTable,
            haystack_normal,
            needle_normal,
        );
    }

    /// Resize pre-allocated buffers to fit a new maximum haystack and
    /// needle size
    pub fn resize(self: *Unicode, max_haystack: usize, max_needle: usize) !void {
        try self.alg.resize(max_haystack, max_needle);
    }

    // Check if buffers have sufficient memory for a given haystack and
    // needle length.
    pub fn hasSize(self: *const Unicode, max_haystack: usize, max_needle: usize) bool {
        return self.alg.hasSize(max_haystack, max_needle);
    }
};

fn doTestScoreUnicode(
    alg: *Unicode,
    haystack: []const u8,
    needle: []const u8,
    comptime score: ?i32,
) !void {
    const s = try alg.score(haystack, needle);
    try std.testing.expectEqual(score, s.?);
}

test "Unicode search" {
    const o = Unicode.Scores{};

    var alg = try Unicode.init(
        std.testing.allocator,
        128,
        32,
        .{},
    );
    defer alg.deinit();

    try doTestScoreUnicode(&alg, "zig⚡ fast", "⚡", o.score_match);
}
