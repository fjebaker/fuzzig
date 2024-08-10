const std = @import("std");

pub const BenchmarkOptions = struct {
    trials: u32 = 10_000,
    warmup: u32 = 100,
};

pub const BenchmarkResult = struct {
    const Self = @This();

    alloc: std.mem.Allocator,
    opts: BenchmarkOptions,
    mean: u64,

    pub fn deinit(_: *Self) void {}

    pub fn printSummary(self: *const Self) void {
        const print = std.debug.print;
        print(
            \\ Benchmark summary for {d} trials:
            \\ Mean: {s}
            \\
        , .{
            self.opts.trials,
            std.fmt.fmtDuration(self.mean),
        });
    }
};

fn invoke(comptime func: anytype, args: std.meta.ArgsTuple(@TypeOf(func))) void {
    const ReturnType = @typeInfo(@TypeOf(func)).Fn.return_type.?;
    switch (@typeInfo(ReturnType)) {
        .ErrorUnion => {
            _ = @call(.never_inline, func, args) catch {
                // std.debug.panic("Benchmarked function returned error {s}", .{err});
            };
        },
        else => _ = @call(.never_inline, func, args),
    }
}

pub fn benchmark(
    alloc: std.mem.Allocator,
    comptime func: anytype,
    args: std.meta.ArgsTuple(@TypeOf(func)),
    opts: BenchmarkOptions,
) !BenchmarkResult {
    var count: usize = 0;
    while (count < opts.warmup) : (count += 1) {
        std.mem.doNotOptimizeAway(true);
        invoke(func, args);
    }
    var timer = try std.time.Timer.start();
    while (count < opts.trials) : (count += 1) {
        std.mem.doNotOptimizeAway(true);
        invoke(func, args);
    }
    const mean = @divFloor(timer.lap(), opts.trials);
    return .{ .alloc = alloc, .opts = opts, .mean = mean };
}
