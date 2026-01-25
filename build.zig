const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const zg = b.lazyDependency("zg", .{
        .target = target,
        .optimize = optimize,
    });

    const mod = b.addModule("fuzzig", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    if (zg) |zg_dep| {
        mod.addImport("code_point", zg_dep.module("code_point"));
        mod.addImport("GeneralCategories", zg_dep.module("GeneralCategories"));
        mod.addImport("LetterCasing", zg_dep.module("LetterCasing"));
        mod.addImport("Normalize", zg_dep.module("Normalize"));
        mod.addImport("CaseFolding", zg_dep.module("CaseFolding"));
    }

    const lib_unit_tests = b.addTest(.{
        .root_module = mod,
    });

    const run_lib_unit_tests = b.addRunArtifact(lib_unit_tests);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_lib_unit_tests.step);

    const exe = b.addExecutable(.{
        .name = "fuzzig-bench",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "fuzzig", .module = mod },
            },
        }),
    });

    const run_cmd = b.addRunArtifact(exe);
    const benchmark_step = b.step("benchmark", "Run benchmarks.");
    benchmark_step.dependOn(&run_cmd.step);
}
